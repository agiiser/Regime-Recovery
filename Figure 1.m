d = 1.0/(sqrt(2*pi));
T = 0.1; % MATURITY
N=51;  
dt=T/(N-1);
u=zeros(51,400,3);
M=400;
Max=1.5;
dx= Max/M;
P= [0,2/3,1/3;0.5,0,0.5;1/3,2/3,0];
lambda=[10,20,10];
mu= [0.08 0.09 0.1];
R=[0.05,0.05,0.05];
SIG=[0.2,0.3,0.4];
st=1;
w=zeros(N);
w(1,1)=1/2;
for i=2:N
    w(i,1)=1/3;
    if mod(i,2)==0
        w(i,i)=4/3;
        w(i-1,i)=1/2;
    else
        w(i,i)=5/6;
        w(i-1,i)=1/3;
    end
    for j=i+1:N
        if mod(i,2)==0
            w(j,i)=4/3;
        else
            w(j,i)=2/3;
        end
    end
end
for m0=2:M
                        if mod(m0,2)==0
                            q(m0)=2/3;
                        else 
                            q(m0)=4/3;
                        end
    end
    q(1)=1/3;
    q(M+1)=1/3;    


tic

    for k=1:3
        for i=1:N
            for j=1:M
                t=(i-1)*dt;
                s=j*dx;
                eta(i,j,k)=s*normcdf((log(s/st)+(R(k)+0.5*SIG(k)^2)*t)/(SIG(k)*sqrt(t)))-st*exp(-R(k)*t)*normcdf((log(s/st)+(R(k)-0.5*SIG(k)^2)*t)/(SIG(k)*sqrt(t)));
            end
        end
    end
    %at maturity price
    for j=1:M
        for i=1:3
            u(1,j,i)=max(0.0,dx*(j)-st);
            eta(1,j,i)=u(1,j,i);
        end
    end
    %Construction of Matrix G( m,m0,l,i)
    for i=1:3
        for l=2:N
            for m=1:M
                for m0=1:M
                    
                    G(m,m0,l,i)=exp(-0.5*(( log(m0/m)-(R(i)-0.5*SIG(i)^2)*((l-1)*dt))/(SIG(i)*sqrt((l-1)*dt) ))^2)*d/(m0*SIG(i)*sqrt((l-1)*dt));
                end
                
            end
        end
    end
    
    
    
     for n=2:N
         
        for m=1:M
            for i=1:3
                vint=0;
                for l=2:n
                    xint=0;
                    
                    for m0=1:M
                        psum=0;
                        for j=1:3
                           psum=psum+P(i,j)*((u(n-l+1,m0,j)-m0*dx));
                        end
                        xint=xint+psum*G(m,m0,l,i)*q(m0+1);
                    end
                    vint=vint+w(n-1,l)*exp(-R(i)*(l-1)*dt)*lambda(i)*exp(-lambda(i)*(l-1)*dt)*(xint+m*dx*exp(R(i)*(l-1)*dt)-st*exp(-R(i)*(n-l)*dt)*(1-logncdf(Max,log(m*dx)+(R(i)-0.5*SIG(i).^2)*(l-1)*dt,SIG(i)*sqrt((l-1)*dt))));
                end
                B(i)=dt*vint+eta(n,m,i)*exp(-lambda(i)*n*dt)-0.5*lambda(i)*dx*dt;
            end
            u(n,m,:)=(eye(3)-w(n-1,1)*dt*diag(lambda)*P)\B';
            u(n,m,:)=max(0.0,u(n,m,:));
        end
    end   
toc



mm=dx:dx:Max;
mmm=1:1:length(mm);
plot(mm,u(N,mmm,1),'b')
hold on
plot(mm,u(N,mmm,2),'g')
hold on
plot(mm,u(N,mmm,3),'r')
plot(mm,eta(N,mmm,1),'--')
hold on
plot(mm,eta(N,mmm,2),'--')
hold on
plot(mm,eta(N,mmm,3),'--')
plot(mm,u(1,mmm,3))
legend('MMGBM State 1','MMGBM State 2','MMGBM State 3','BSM State 1','BSM State 2','BSM State 3','Maturity Value')