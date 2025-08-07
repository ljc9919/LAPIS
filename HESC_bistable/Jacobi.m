function [Jacob]=Jacobi(mu,IM,C,D,S,Hilln,G,g,k)
kk = size(mu,2);
H0 = -diag(k')*eye(kk);
H=zeros(kk,kk);
for i=1:kk
    for j=1:kk
        if IM(i,j)==1
            H(j,i)=(G(j)-g(j))*dHa(mu(i),S(i,j),C(i,j),Hilln(i,j));
        end
        if IM(i,j)==-1
            H(j,i)=-g(j)*dHa(mu(i),S(i,j),D(i,j),Hilln(i,j));
        end
    end
end
Jacob=H0+H;
end