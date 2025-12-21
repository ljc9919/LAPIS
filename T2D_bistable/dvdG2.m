function [dv1dG,dv2dG] = dvdG2(g,G,k,mu1,mu2,S,C,D,N)
A = zeros(N,N);
i = 1;
    A(i,i) = -k(i)-g(i)*(C(i,i)+D(i,i))*4*mu1(i)^3*S(i,i)^4/(S(i,i)^4+mu1(i)^4)^2+...
             G(i)*C(i,i)*4*mu1(i)^3*S(i,i)^4/(S(i,i)^4+mu1(i)^4)^2;          
    for j = setdiff(1:N,i)
        A(i,j) = -g(i)*(C(j,i)+D(j,i))*4*mu1(j)^3*S(j,i)^4/(S(j,i)^4+mu1(j)^4)^2+...
                 G(i)*C(j,i)*4*mu1(j)^3*S(j,i)^4/(S(j,i)^4+mu1(j)^4)^2;  
    end
    for m = setdiff(1:N,i)
        A(m,m) = -k(m)-g(m)*(C(m,m)+D(m,m))*4*mu1(m)^3*S(m,m)^4/(S(m,m)^4+mu1(m)^4)^2+...
                 G(m)*C(m,m)*4*mu1(m)^3*S(m,m)^4/(S(m,m)^4+mu1(m)^4)^2; 
             for q = setdiff(1:N,m)
                 A(m,q) = -g(m)*(C(q,m)+D(q,m))*4*mu1(q)^3*S(q,m)^4/(S(q,m)^4+mu1(q)^4)^2+...
                 G(m)*C(q,m)*4*mu1(q)^3*S(q,m)^4/(S(q,m)^4+mu1(q)^4)^2; 
             end
    end
    
    for ii = 1:5
        b = zeros(N,1);
        b(ii) = -sum(C(:,ii)'.*mu1.^4./(S(:,ii)'.^4+mu1.^4));  
        dv1dG(:,ii) = A\b;
    end
 
A = zeros(N,N);
i = 1;
    A(i,i) = -k(i)-g(i)*(C(i,i)+D(i,i))*4*mu2(i)^3*S(i,i)^4/(S(i,i)^4+mu2(i)^4)^2+...
             G(i)*C(i,i)*4*mu2(i)^3*S(i,i)^4/(S(i,i)^4+mu2(i)^4)^2;          
    for j = setdiff(1:N,i)
        A(i,j) = -g(i)*(C(j,i)+D(j,i))*4*mu2(j)^3*S(j,i)^4/(S(j,i)^4+mu2(j)^4)^2+...
                 G(i)*C(j,i)*4*mu2(j)^3*S(j,i)^4/(S(j,i)^4+mu2(j)^4)^2;  
    end
    for m = setdiff(1:N,i)
        A(m,m) = -k(m)-g(m)*(C(m,m)+D(m,m))*4*mu2(m)^3*S(m,m)^4/(S(m,m)^4+mu2(m)^4)^2+...
                 G(m)*C(m,m)*4*mu2(m)^3*S(m,m)^4/(S(m,m)^4+mu2(m)^4)^2; 
             for q = setdiff(1:N,m)
                 A(m,q) = -g(m)*(C(q,m)+D(q,m))*4*mu2(q)^3*S(q,m)^4/(S(q,m)^4+mu2(q)^4)^2+...
                 G(m)*C(q,m)*4*mu2(q)^3*S(q,m)^4/(S(q,m)^4+mu2(q)^4)^2; 
             end
    end

    for ii = 1:5
        b = zeros(N,1);
        b(ii) = -sum(C(:,ii)'.*mu2.^4./(S(:,ii)'.^4+mu2.^4));  
        dv2dG(:,ii) = A\b;
    end
end
