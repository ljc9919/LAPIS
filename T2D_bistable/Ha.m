function H=Ha(X,S,a,n)
H=a.*X.^n./(S.^n+X.^n);
end