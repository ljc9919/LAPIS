function [dv1ds,dv2ds,dv3ds] = dvds(g,G,k,mu1,mu2,mu3,S,C,D,N)
A = zeros(N,N);
for i = 1:N
    b = zeros(N,1);
    A(i,i) = -k(i)-g(i)*(C(i,i)+D(i,i))*4*mu1(i)^3*S(i,i)^4/(S(i,i)^4+mu1(i)^4)^2+...
             G(i)*C(i,i)*4*mu1(i)^3*S(i,i)^4/(S(i,i)^4+mu1(i)^4)^2;          
    b(i) = ((G(i)-g(i))*C(i,i)*4*mu1(i)^4*S(i,i)^3)/((S(i,i)^4+mu1(i)^4)^2)-...
           (g(i)*D(i,i)*4*mu1(i)^4*S(i,i)^3)/((S(i,i)^4+mu1(i)^4)^2);
    for j = setdiff(1:N,i)
        A(i,j) = -g(i)*(C(j,i)+D(j,i))*4*mu1(j)^3*S(j,i)^4/(S(j,i)^4+mu1(j)^4)^2+...
                 G(i)*C(j,i)*4*mu1(j)^3*S(j,i)^4/(S(j,i)^4+mu1(j)^4)^2;  
    end
    for m = setdiff(1:N,i)
        A(m,m) = -k(m)-g(m)*(C(m,m)+D(m,m))*4*mu1(m)^3*S(m,m)^4/(S(m,m)^4+mu1(m)^4)^2+...
                 G(m)*C(m,m)*4*mu1(m)^3*S(m,m)^4/(S(m,m)^4+mu1(m)^4)^2; 
        b(m) = ((G(m)-g(m))*C(i,m)*4*mu1(i)^4*S(i,m)^3)/((S(i,m)^4+mu1(i)^4)^2)-...
               (g(m)*D(i,m)*4*mu1(i)^4*S(i,m)^3)/((S(i,m)^4+mu1(i)^4)^2);   
             for q = setdiff(1:N,m)
                 A(m,q) = -g(m)*(C(q,m)+D(q,m))*4*mu1(q)^3*S(q,m)^4/(S(q,m)^4+mu1(q)^4)^2+...
                          G(m)*C(q,m)*4*mu1(q)^3*S(q,m)^4/(S(q,m)^4+mu1(q)^4)^2; 
             end
    end
    dv1ds(:,i) = A\b;
end

A = zeros(N,N);
for i = 1:N
    b = zeros(N,1);
    A(i,i) = -k(i)-g(i)*(C(i,i)+D(i,i))*4*mu2(i)^3*S(i,i)^4/(S(i,i)^4+mu2(i)^4)^2+...
             G(i)*C(i,i)*4*mu2(i)^3*S(i,i)^4/(S(i,i)^4+mu2(i)^4)^2;          
    b(i) = ((G(i)-g(i))*C(i,i)*4*mu2(i)^4*S(i,i)^3)/((S(i,i)^4+mu2(i)^4)^2)-...
           (g(i)*D(i,i)*4*mu2(i)^4*S(i,i)^3)/((S(i,i)^4+mu2(i)^4)^2);
    for j = setdiff(1:N,i)
        A(i,j) = -g(i)*(C(j,i)+D(j,i))*4*mu2(j)^3*S(j,i)^4/(S(j,i)^4+mu2(j)^4)^2+...
                 G(i)*C(j,i)*4*mu2(j)^3*S(j,i)^4/(S(j,i)^4+mu2(j)^4)^2;  
    end
    for m = setdiff(1:N,i)
        A(m,m) = -k(m)-g(m)*(C(m,m)+D(m,m))*4*mu2(m)^3*S(m,m)^4/(S(m,m)^4+mu2(m)^4)^2+...
                 G(m)*C(m,m)*4*mu2(m)^3*S(m,m)^4/(S(m,m)^4+mu2(m)^4)^2; 
        b(m) = ((G(m)-g(m))*C(i,m)*4*mu2(i)^4*S(i,m)^3)/((S(i,m)^4+mu2(i)^4)^2)-...
               (g(m)*D(i,m)*4*mu2(i)^4*S(i,m)^3)/((S(i,m)^4+mu2(i)^4)^2);   
             for q = setdiff(1:N,m)
                 A(m,q) = -g(m)*(C(q,m)+D(q,m))*4*mu2(q)^3*S(q,m)^4/(S(q,m)^4+mu2(q)^4)^2+...
                          G(m)*C(q,m)*4*mu2(q)^3*S(q,m)^4/(S(q,m)^4+mu2(q)^4)^2; 
             end
    end
    dv2ds(:,i) = A\b;
end

A = zeros(N,N);
for i = 1:N
    b = zeros(N,1);
    A(i,i) = -k(i)-g(i)*(C(i,i)+D(i,i))*4*mu3(i)^3*S(i,i)^4/(S(i,i)^4+mu3(i)^4)^2+...
             G(i)*C(i,i)*4*mu3(i)^3*S(i,i)^4/(S(i,i)^4+mu3(i)^4)^2;          
    b(i) = ((G(i)-g(i))*C(i,i)*4*mu3(i)^4*S(i,i)^3)/((S(i,i)^4+mu3(i)^4)^2)-...
           (g(i)*D(i,i)*4*mu3(i)^4*S(i,i)^3)/((S(i,i)^4+mu3(i)^4)^2);
    for j = setdiff(1:N,i)
        A(i,j) = -g(i)*(C(j,i)+D(j,i))*4*mu3(j)^3*S(j,i)^4/(S(j,i)^4+mu3(j)^4)^2+...
                 G(i)*C(j,i)*4*mu3(j)^3*S(j,i)^4/(S(j,i)^4+mu3(j)^4)^2;  
    end
    for m = setdiff(1:N,i)
        A(m,m) = -k(m)-g(m)*(C(m,m)+D(m,m))*4*mu3(m)^3*S(m,m)^4/(S(m,m)^4+mu3(m)^4)^2+...
                 G(m)*C(m,m)*4*mu3(m)^3*S(m,m)^4/(S(m,m)^4+mu3(m)^4)^2; 
        b(m) = ((G(m)-g(m))*C(i,m)*4*mu3(i)^4*S(i,m)^3)/((S(i,m)^4+mu3(i)^4)^2)-...
               (g(m)*D(i,m)*4*mu3(i)^4*S(i,m)^3)/((S(i,m)^4+mu3(i)^4)^2);   
             for q = setdiff(1:N,m)
                 A(m,q) = -g(m)*(C(q,m)+D(q,m))*4*mu3(q)^3*S(q,m)^4/(S(q,m)^4+mu3(q)^4)^2+...
                          G(m)*C(q,m)*4*mu3(q)^3*S(q,m)^4/(S(q,m)^4+mu3(q)^4)^2; 
             end
    end
    dv3ds(:,i) = A\b;
end


end
