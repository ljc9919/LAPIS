function dG2dfdx = dG2_dfdx(g,G,mu,S,C,D,N,dv1dG)
dG2dfdx = zeros(N,N,N);
for i = 1:N
    for j = 1:N
        for p = 1:N
            if p == i
                dG2dfdx(p,j,i) = C(j,p)*4*mu(j)^3*S(j,p)^4/(mu(j)^4+S(j,p)^4)^2+...
  G(p)*C(j,p)*4*S(j,p)^4*(3*mu(j)^2*dv1dG(j,p)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dG(j,p))/(mu(j)^4+S(j,p)^4)^4-...
  g(p)*(C(j,p)+D(j,p))*4*S(j,p)^4*(3*mu(j)^2*dv1dG(j,p)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dG(j,p))/(mu(j)^4+S(j,p)^4)^4;
            elseif p~=i
                dG2dfdx(p,j,i) = ...
  G(p)*C(j,p)*4*S(j,p)^4*(3*mu(j)^2*dv1dG(j,i)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dG(j,i))/(mu(j)^4+S(j,p)^4)^4-...
  g(p)*(C(j,p)+D(j,p))*4*S(j,p)^4*(3*mu(j)^2*dv1dG(j,i)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dG(j,i))/(mu(j)^4+S(j,p)^4)^4;
            end            
        end
    end
end