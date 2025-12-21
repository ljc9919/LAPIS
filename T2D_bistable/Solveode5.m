function [mu,sigma,alpha] = Solveode5(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax)
[xx,n,sigma]=Solver(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax);
N = 5;
index=size(n,1);  %% The number of the stable states
alpha=zeros(index,1);  %% The weight of the stable states
mu=zeros(index,N);  %% The mean value of the Gaussian density function
for i=1:index
%The mean value of each stable state
    mu(i,:)=xx(n(i,1),:); 
    %sigma0{i}=reshape(sigma(:,:,i),N,N)';
%The weight of each stable state
    alpha(i)=n(i,2)/sum(n(:,2));   %%
end
end