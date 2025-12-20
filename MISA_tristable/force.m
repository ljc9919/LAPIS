function f=force(t,x,k,G,g,S,IM,Hilln,AB)
[num,m]=size(x);
H1=zeros(num,m,num); 
H2=zeros(num,m,num); 
for i=1:2
    for j=1:2
        if IM(i,j)==1
            H1(i,:,j)=Ha(x(i,:),S(i,j),AB(i,j),Hilln(i,j));%activation
        elseif IM(i,j)==-1
            H2(i,:,j)=Ha(x(i,:),S(i,j),AB(i,j),Hilln(i,j));%inhibition
        end
    end
end
for i=1:2
    f(i,:)=g(i)+G(i)*sum(H1(:,:,i))-g(i)*sum(H1(:,:,i))-g(i)*sum(H2(:,:,i))-k(i)*x(i,:);
end
end
function H=Ha(X,S,a,n)
H=a.*X.^n./(S.^n+X.^n);
end
   

