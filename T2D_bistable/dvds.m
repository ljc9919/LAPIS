function [dv1ds,dv2ds] = dvds(g,G,k,mu1,mu2,S,C,D,N)
A = zeros(N,N);
i = 1;
    A(i,i) = -k(i)-g(i)*(C(i,i)+D(i,i))*4*mu1(i)^3*S(i,i)^4/(S(i,i)^4+mu1(i)^4)^2+...
             G(i)*C(i,i)*4*mu1(i)^3*S(i,i)^4/(S(i,i)^4+mu1(i)^4)^2;          
    for j = setdiff(1:N,i)
        A(i,j) = -g(i)*(C(j,i)+D(j,i))*4*mu1(j)^3*S(j,i)^4/(S(j,i)^4+mu1(j)^4)^2+...
                 G(i)*C(j,i)*4*mu1(j)^3*S(j,i)^4/(S(j,i)^4+mu1(j)^4)^2;  
    end
    for m = setdiff(1:N,i)
        A(m,m) = -k(m)-g(m)*(C(m,m)+D(m,m))*4*mu1(m)^3*S(m,m)^4/(S(m,m)^4+mu1(m)^4)^2+...    %C(i,m)+D(i,m)
                 G(m)*C(m,m)*4*mu1(m)^3*S(m,m)^4/(S(m,m)^4+mu1(m)^4)^2; 
             for q = setdiff(1:N,m)
                 A(m,q) = -g(m)*(C(q,m)+D(q,m))*4*mu1(q)^3*S(q,m)^4/(S(q,m)^4+mu1(q)^4)^2+...
                          G(m)*C(q,m)*4*mu1(q)^3*S(q,m)^4/(S(q,m)^4+mu1(q)^4)^2; 
             end
    end
    
    for ii = 1:5
        b = zeros(N,1);
        b(ii) = ((G(ii)-g(ii))*C(ii,ii)*4.*mu1(ii).^4.*S(ii,ii)'.^3)/((S(ii,ii)'.^4+mu1(ii).^4).^2)-...
         g(ii)*D(ii,ii).*4.*mu1(ii).^4.*S(ii,ii)'.^3./((S(ii,ii)'.^4+mu1(ii).^4).^2);
        for jj = setdiff(1:N,i)
            b(jj) = (G(jj)-g(jj))*C(ii,jj)*4.*mu1(ii).^4.*S(ii,jj)'.^3./((S(ii,jj)'.^4+mu1(ii).^4).^2)-...
         g(jj)*D(ii,jj)'.*4.*mu1(ii).^4.*S(ii,jj)'.^3./((S(ii,jj)'.^4+mu1(ii).^4).^2);
        end
        dv1ds(:,ii) = A\b;
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
        b(ii) = ((G(ii)-g(ii))*C(ii,ii)*4.*mu2(ii).^4.*S(ii,ii)'.^3)/((S(ii,ii)'.^4+mu2(ii).^4).^2)-...
         g(ii)*D(ii,ii).*4.*mu2(ii).^4.*S(ii,ii)'.^3./((S(ii,ii)'.^4+mu2(ii).^4).^2);
        for jj = setdiff(1:N,i)
            b(jj) = (G(jj)-g(jj))*C(ii,jj)*4.*mu2(ii).^4.*S(ii,jj)'.^3./((S(ii,jj)'.^4+mu2(ii).^4).^2)-...
         g(jj)*D(ii,jj)'.*4.*mu2(ii).^4.*S(ii,jj)'.^3./((S(ii,jj)'.^4+mu2(ii).^4).^2);
        end
        dv2ds(:,ii) = A\b;
    end
    
end
