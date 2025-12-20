function [Sig]=Sigma(xx,IM,C,D,S,Hilln,G,g,K,d)
% xx is the vector of sigma which is Row Major Order
%IM: Interaction matrix.
%AB:the sthength of activation and inhibution.
%S:Threshold matrix.
%n:Hill function matrix.

kk=length(xx); %the dimension of the system.

H0=-diag(K')*eye(kk);
H=zeros(kk,kk);
for i=1:kk
    for j=1:kk
        if IM(i,j)==1
            H(j,i)=(G(j)-g(j))/max(sum(C))*Ha(xx(i),S(i,j),C(i,j),Hilln(i,j));
        end
        if IM(i,j)==-1
            H(j,i)=-g(j)/max(sum(D))*Ha(xx(i),S(i,j),D(i,j),Hilln(i,j));
        end
    end
end
H=H0+H;
% A*sigma+sigma*A'+2D=0
% a=1;b=0.5;S=0.7;n=4;k=1;
% syms A B C; 
% Ajac=jacobian([Ha0(A,S,a,n)+Hr0(B,S,b,n)-k*A;
% Hr0(C,S,b,n)+Hr0(A,S,b,n)-k*B;
% Ha0(A,S,a,n)-k*C],[A,B,C]);
% Ajac=subs(Ajac,{'A','B','C'},{xx(1),xx(2),xx(3)});
% Ajac=double(Ajac);

P=zeros(kk^2,kk^2);  %coefficient matrix

%%the initial of coeffiicient matrix
for i=0:(kk-1)
    P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)=P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)+H;
end

for m=0:kk-1
    for i=1:kk
        for j=1:kk
            P(m*kk+i,(j-1)*kk+i)=P(m*kk+i,(j-1)*kk+i)+H(m+1,j);
        end
    end
end

B=zeros(kk^2,1);
for i=1:kk
    B((i-1)*kk+i)=-2*d;
end

Sig=P\B;
Sig=reshape(Sig,kk,kk)';
end
function H=Hr(X,S,b,n)
H=-b*n*S^n*X^(n-1)./(X^n+S^n)^2;
end
function H=Ha(X,S,a,n)
H=a*n*S^n*X^(n-1)./(S^n+X^n)^2;
end
