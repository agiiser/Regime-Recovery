clear all
IV=xlsread('estimation.xlsx','D2:D1403');
%IV=xlsread('IV4.xlsx','D2:D201');
A=length(IV);
for i=2:A
    delta(i)=IV(i)-IV(i-1);
end
B=2*prctile(abs(delta),90);
flat(1)=IV(1);
M(1)=0;
N(1)=IV(1);
P(1)=1;
R(1)=IV(1);
for i=2:A
    if abs(delta(i))<B
        flat(i)=flat(i-1);
    else
        flat(i)=flat(i-1)+delta(i);
    end
    if mod(i-1,20)==0
        M(i)=M(i-1)+flat(i)-flat(i-1);
    else
        M(i)=M(i-1);
    end
    N(i)=flat(i)-M(i);
    if N(i-1)==N(i)
        P(i)=P(i-1);
    else
        P(i)=i;
    end
    if P(i)==P(i-1)
        R(i)=min(R(i-1),IV(i));
    else
        R(i)=IV(P(i));
    end
    if P(i)==P(i-1)
        R2(i)=R(i-1)-R(i);
    else
        R2(i)=0;
    end
end

a=2*prctile(R2,98);
R3(A)=A;
for i=1:A-1
    if N(i+1)~=N(i)
        R3(i)=i;
    else
        R3(i)=0;
    end
       
end

for j=1:A
    for i=1:A
        if R3(i)==0
           if R2(i)<a
              R3(i)=R3(i+1);
           else
              R3(i)=i;
           end
        end
    end
end


for i=2:A
    if abs(R2(i))<a
        b=1000;
        for k=i:R3(i)
            b=min(b,R(k));
        end
        S3(i)=b;
    else
        S3(i)=R(i);
    end
end
S3(1)=S3(2);
c=(max(S3)-min(S3))/10;
d=min(S3);
e=max(S3);
f=prctile(S3,50);
for i=2:A
    if abs(S3(i)-d)<c
        T(i)=d;
    elseif abs(S3(i)-e)<c
        T(i)=e;
    else
        T(i)=f;
    end
end
hist(IV,10)








        

    


