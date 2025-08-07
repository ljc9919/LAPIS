function f=force2(t,x,k,G,g,S,IM,Hilln,AB,C,D)
%This is an implementation of solving the stable states of the ODEs
%%x:current location
%%a, aa, b, sa, sb, k, n: the parameters
%%mat:the matrix representing the relationship between genes
[num,m]=size(x);
H1=zeros(num,m,num); 
H2=zeros(num,m,num); 
for i=1:52
    for j=1:52
        if IM(i,j)==1
            H1(i,:,j)=Ha(x(i,:),S(i,j),AB(i,j),Hilln(i,j));
        elseif IM(i,j)==-1
            H2(i,:,j)=Ha(x(i,:),S(i,j),AB(i,j),Hilln(i,j));
        end
    end
end
for i=1:52
    f(i,:)=g(i)+G(i)/max(sum(C))*sum(H1(:,:,i))-g(i)/max(sum(C))*sum(H1(:,:,i))-g(i)/max(sum(D))*sum(H2(:,:,i))-k(i)*x(i,:);
end
end
function H=Hr(X,S,b,n)
H=b.*S.^n./(X.^n+S.^n);
end
function H=Ha(X,S,a,n)
H=a.*X.^n./(S.^n+X.^n);
end
   

