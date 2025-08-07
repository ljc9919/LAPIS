clear,clc
tic()
IM=load ('matrixCF.txt');
N = size(IM,1); %The dimension of the system.
Hilln=4*ones(N,N);

AB1 = IM; AB = AB1; AB(find(AB1<0)) = -AB(find(AB1<0));

cno = (AB1(:)>0);
C = zeros(N,N);
C(cno) = AB1(cno);

dno = (AB1(:)<0);
D = zeros(N,N);
D(dno) = -AB1(dno);

k=ones(N,1);
wemax = 5*ones(1,N);
cycle_index = 1000;
dd = 0.01;

G = zeros(N,1);
g = zeros(N,1);
S = zeros(N,N);
no = 1;

        G(1:11)=15;
        G(12:22)=9;
        G(23:end)=5;
        g(1:11)=1;
        g(12:22)=9;
        g(23:end)=3;
        S(1:11,:)=4;
        S(12:22,:)=2;
        S(23:end,:)=3;

k = 1*ones(N,1);

[xx,n,sigma]=Solver(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax);

N = 52;
index=size(n,1);  %% The number of the stable states
alpha=zeros(index,1);  %% The weight of the stable states
mu=zeros(index,N);  %% The mean value of the Gaussian density function
for i=1:index
%The mean value of each stable state
    mu(i,:)=xx(n(i,1),:);
    sigma0{i}=reshape(sigma(:,:,i),N,N)';
%The weight of each stable state
    alpha(i)=n(i,2)/sum(n(:,2));   %%%%%????a??alpha?¡§???¡§¡§¡ê¡è?¨º?muo¡§aalpha¡§o???¡§?||¨¬???¨º
end
     
   

save a52searchG2is1G3is0.5g3is0.1S2andS30.5_.mat