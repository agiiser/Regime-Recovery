d = 1.0/(sqrt(2*pi));
T = 0.2; % MATURITY
N=51;
dt=T/(N-1);

M=400;
Max=3;
dx= Max/M;
P= [0,2/3,1/3;0.5,0,0.5;1/3,2/3,0];
lambda=[10,20,10];
mu= [0.08 0.09 0.1];
R=[0,0,0];
SIG=[0.2,0.3,0.4];
X(1)=1;
S(1)=1;
tic
for p=2:200
    p1=1-lambda(X(p-1))*dt;
    k2=rand;
    if k2<=p1
        X(p)=X(p-1);
    elseif X(p-1)==1
        p11=P(1,2);
        k11=rand;
        if k11<=p11
            X(p)=2;
        else
            X(p)=3;
        end
    elseif X(p-1)==2
        p11=P(2,1);
        k11=rand;
        if k11<=p11
            X(p)=1;
        else
            X(p)=3;
        end
    else
        p11=P(3,2);
        k11=rand;
        if k11<=p11
            X(p)=2;
        else
            X(p)=1;
        end
    end
    %simulation of stock prices
    if X(p)==1
        S(p)=S(p-1)*exp((mu(1)-0.5*SIG(1)^2)*dt+SIG(1)*sqrt(dt)*randn);
    elseif X(p)==2
        S(p)=S(p-1)*exp((mu(2)-0.5*SIG(2)^2)*dt+SIG(2)*sqrt(dt)*randn);
    else
        S(p)=S(p-1)*exp((mu(3)-0.5*SIG(3)^2)*dt+SIG(3)*sqrt(dt)*randn);
    end
end
%decreasing time to maturity
for p=1:200
    L=[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280];
    for i=2:14
        if (p<=L(i))&&(p>L(i-1))
            Y(p)=L(i+1)-p;
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



for p=1:200
    
    
    st=round(S(p),2);
    %BSM option price
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
    
       
    %option prices at grid points  
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
                B(i)=dt*vint+eta(n,m,i)*exp(-lambda(i)*n*dt);
            end
            u(n,m,:)=(eye(3)-w(n-1,1)*dt*diag(lambda)*P)\B';
            u(n,m,:)=max(0.0,u(n,m,:));
        end
    end   
    %Interpolation
    A=S(p)/dx;
    
    if X(p)==1
        C(p)=((A-floor(A))*u(Y(p),ceil(A),1))+((floor(A+1)-A)*u(Y(p),floor(A),1));
        
    elseif X(p)==2
        C(p)=((A-floor(A))*u(Y(p),ceil(A),2))+((floor(A+1)-A)*u(Y(p),floor(A),2));
        
    else
        C(p)=((A-floor(A))*u(Y(p),ceil(A),3))+((floor(A+1)-A)*u(Y(p),floor(A),3));
        
    end
    %Implied volatility
    IV(p)=blsimpv(S(p),st,0,(Y(p)*dt),C(p));
    D(p)=lognrnd(-.04,0.3)*IV(p);
end
toc
xlswrite('IV4.xlsx',transpose(S),'A2:A201')
xlswrite('IV4.xlsx',transpose(C),'B2:B201')
xlswrite('IV4.xlsx',transpose(X),'C2:C201')
xlswrite('IV4.xlsx',transpose(IV),'D2:D201')
xlswrite('IV4.xlsx',transpose(SIG(X)),'E2:E201')
plot(IV,'.')
hold on
plot(SIG(X),'*')


    






    
        
