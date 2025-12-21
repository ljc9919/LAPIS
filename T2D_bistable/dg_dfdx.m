function dgdfdx = dg_dfdx(g,G,mu,S,C,D,N,dv1dg)
dgdfdx = zeros(N,N,N);
for i = 1:N
    for j = 1:N
        for p = 1:N
            if  p == i
                dgdfdx(p,j,i) = -(C(j,p)+D(j,p))*4*mu(j)^3*S(j,p)^4/(mu(j)^4+S(j,p)^4)^2+...
  G(p)*C(j,p)*4*S(j,p)^4*(3*mu(j)^2*dv1dg(j,p)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dg(j,p))/(mu(j)^4+S(j,p)^4)^4-...
  g(p)*(C(j,p)+D(j,p))*4*S(j,p)^4*(3*mu(j)^2*dv1dg(j,p)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dg(j,p))/(mu(j)^4+S(j,p)^4)^4;
            else
                dgdfdx(p,j,i) = ...
  G(p)*C(j,p)*4*S(j,p)^4*(3*mu(j)^2*dv1dg(j,i)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dg(j,i))/(mu(j)^4+S(j,p)^4)^4-...
  g(p)*(C(j,p)+D(j,p))*4*S(j,p)^4*(3*mu(j)^2*dv1dg(j,i)*(mu(j)^4+S(j,p)^4)^2-mu(j)^6*8*(mu(j)^4+S(j,p)^4)*dv1dg(j,i))/(mu(j)^4+S(j,p)^4)^4;
            end            
        end
    end
end
  
     
                            
                             
                          
