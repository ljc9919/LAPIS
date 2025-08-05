function dsdfdx = ds_dfdx(g,G,mu,S,C,D,N,dv1ds)
dsdfdx = zeros(N,N,N);
for i = 1:N
    for j = 1:N
        for p = 1:N
            if j == i
                dsdfdx(p,j,i) = ...
  (G(p)*C(j,p)-g(p)*(C(j,p)+D(j,p)))*4*(-5*mu(j)^10*dv1ds(j,j)*S(j,p)^4+3*mu(j)^2*dv1ds(j,j)*S(j,p)^12-2*mu(j)^6*dv1ds(j,j)*S(j,p)^8+4*mu(j)^11*S(j,p)^3-4*mu(j)^3*S(j,p)^11)/(mu(j)^4+S(j,p)^4)^4;          
            elseif j~=i
                dsdfdx(p,j,i) = ...
  (G(p)*C(j,p)-g(p)*(C(j,p)+D(j,p)))*4*S(j,p)^4*(-5*mu(j)^10*dv1ds(j,i)+3*mu(j)^2*dv1ds(j,i)*S(j,p)^8-2*mu(j)^6*dv1ds(j,i)*S(j,p)^4)/(mu(j)^4+S(j,p)^4)^4;
            end            
        end
    end
end