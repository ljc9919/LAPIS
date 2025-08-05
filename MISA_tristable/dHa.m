function H=dHa(X,S,a,n)
H=a*n*S^n*X^(n-1)./(S^n+X^n)^2;
end