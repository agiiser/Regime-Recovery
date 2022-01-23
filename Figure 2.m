d = 1.0/(sqrt(2*pi));
T = 0.2; % MATURITY
N=16;  
dt=T/(N-1);
RR=0.01;
M=400;
Max=3;
dx= Max/M;
P= [0,2/3,1/3;0.5,0,0.5;1/3,2/3,0];
lambda=[10,20,10];
mu= [0.08 0.09 0.1];
R=[RR,RR,RR];
SIG=[0.2,0.3,0.4];
L=floor(M/Max);
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
for i=1:3
        for l=2:N
            for m=1:M
                for m0=1:M
                    
                    G(m,m0,l,i)=exp(-0.5*(( log(m0/m)-(R(i)-0.5*SIG(i)^2)*((l-1)*dt))/(SIG(i)*sqrt((l-1)*dt) ))^2)*d/(m0*SIG(i)*sqrt((l-1)*dt));
                end
                
            end
        end
end
for a=1:21
         st=0.8+0.02*(a-1);
         ST(a)=st;
         if st <= 1
              for k=1:3
                  for i=2:N
                      for j=1:M
                          t=(i-1)*dt;
                          s=(j)*dx;
                          etap(i,j,k)=-s*normcdf(-(log(s/st)+(R(k)+0.5*SIG(k)^2)*t)/(SIG(k)*sqrt(t)))+st*exp(-R(k)*t)*normcdf(-(log(s/st)+(R(k)-0.5*SIG(k)^2)*t)/(SIG(k)*sqrt(t)));
                      end
                  end
              end
              for j=1:M
                  for i=1:3
                      uu(1,j,i)=max(0.0,st-dx*(j-1));
                      etap(1,j,i)=uu(1,j,i);
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
                                      psum=psum+P(i,j)*uu(n-l+1,m0,j);
                                  end
                                  xint=xint+psum*G(m,m0,l,i)*q(m0+1);
                              end
                              vint=vint+ w(n-1,l)*exp(-R(i)*(l-1)*dt)*xint*lambda(i)*exp(-lambda(i)*(l-1)*dt);
                          end
                          B(i)=dt*vint+etap(n,m,i)*exp(-lambda(i)*n*dt);
                      end
                      uu(n,m,:)=(eye(3)-w(n-1,1)*dt*diag(lambda)*P)\B';
                      uu(n,m,:)=max(0.0,uu(n,m,:));
                  end
              end
              for i=1:3
                  Price(a,i)=uu(N,L,i);
                  IV(a,i) = blsimpv(1, ST(a), RR, N*dt, uu(N,L,i),[],[],[],false);
              end
         else
             for k=1:3
                 for i=1:N
                     for j=1:M
                         t=(i-1)*dt;
                         s=j*dx;
                         eta(i,j,k)=s*normcdf((log(s/st)+(R(k)+0.5*SIG(k)^2)*t)/(SIG(k)*sqrt(t)))-st*exp(-R(k)*t)*normcdf((log(s/st)+(R(k)-0.5*SIG(k)^2)*t)/(SIG(k)*sqrt(t)));
                     end
                 end
             end
    
             for j=1:M
                 for i=1:3
                     u(1,j,i)=max(0.0,dx*(j)-st);
                     eta(1,j,i)=u(1,j,i);
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
                         B(i)=dt*vint+eta(n,m,i)*exp(-lambda(i)*n*dt);
                     end
                     u(n,m,:)=(eye(3)-w(n-1,1)*dt*diag(lambda)*P)\B';
                     u(n,m,:)=max(0.0,u(n,m,:));
                 end
             end 
             for i=1:3
                 Price(a,i)=u(N,L,i);
                 IV(a,i) = blsimpv(1, ST(a), RR, N*dt, u(N,L,i),[],[],[],true);
             end
         end
         clear etap;
         clear eta;
         clear g;
end
plot(ST,IV(:,1)','*')
hold on
plot(ST,IV(:,2)','o')
hold on
plot(ST,IV(:,3)','*')
legend('State 1','State 2','State 3')

             
         
                
         
                  