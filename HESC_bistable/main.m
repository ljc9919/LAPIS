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
wemax = 10*ones(1,N);
cycle_index = 1000;
dd = 0.01;

G = zeros(N,1);
g = zeros(N,1);
S = zeros(N,N);
no = 1;

for ss=5:2:15
    for qq=1:2:10
        for kk=1:2:ss
            for ll=1:2:qq
        G(1:11)=ss;
        G(12:22)=qq;
        G(23:end)=3;
        g(1:11)=kk;
        g(12:22)=ll;
        g(23:end)=1;
        S(1:11,:)=3;
        S(12:22,:)=3;
        S(23:end,:)=3;

k = 1*ones(N,1);

[xx,n]=Solver(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax);

N = 52;
index=size(n,1);  %% The number of the stable states
alpha=zeros(index,1);  %% The weight of the stable states
mu=zeros(index,N);  %% The mean value of the Gaussian density function
for i=1:index
%The mean value of each stable state
    mu(i,:)=xx(n(i,1),:); 
%The weight of each stable state
    alpha(i)=n(i,2)/sum(n(:,2));   %%%%%????a??alpha?¡§???¡§¡§¡ê¡è?¨º?muo¡§aalpha¡§o???¡§?||¨¬???¨º
end
     mu(find(alpha<=0.01),:)=[];
     alpha(find(alpha<=0.01))=[];
     if size(mu,1)==2
        loveG{no} = G ;loveg{no}=g;lovemu{no}=mu;
        loveS{no} = S; no=no+1;
     end
    end
    end
        end
    end


save a52searchG15to15G21to10G3is3g3is1Sall13.mat