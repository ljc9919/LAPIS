clear,clc
tic()
mutrue1 = [16.9089	0.0151;11.7317	12.7648;0.0637	14.8343]; %(highlow,highhigh,lowhigh)ground_truth mu value
mutrue = reshape(mutrue1,1,6);
datamu1 = mutrue1(1,:);   datamu2 = mutrue1(2,:);   datamu3 = mutrue1(3,:); 

dataweight = [0.35500 0.38567 0.25933];%ground_truth weight value
datasigma(:,:,1) = [0.4000	-0.0036;-0.0036	0.4184];
datasigma(:,:,2) = [0.4585	-0.0141;-0.0141	0.4406];
datasigma(:,:,3) = [0.4087	-0.0007;-0.0007	0.4000];%ground_truth sigma value

N = 2;%Dim
IM = [1 -1;-1 1]; %prior network structure
Hilln = 4*ones(N,N); %Hill coefficient
AB = IM; AB(find(IM<0)) = -AB(find(IM<0));
cno = (IM(:)>0); C = zeros(N,N); C(cno) = IM(cno);%activation structure
dno = (IM(:)<0); D = zeros(N,N); D(dno) = -IM(dno);%inhibition structure

k = 1*ones(N,1);%fixed degradation value
wemax = 10*ones(1,N);%initial values for ODE drawn from U[0,wemax]
lambda_penalty = 100;
cycle_index = 3000;%sampled initial points for ODE
dd = 0.1; %assigned arbitrarily
min_loss = inf;

for ii = 1:5
    for jj = 1:10
        for ss = 1:5
            for sss = 1:5
                G = [10+ii*2;10+ii*2];%two genes assumed to be identical 12 14 16 18 20
                g = [jj;jj];%two genes assumed to be identical 1 2 3...10
                S = [ss*2,ss*2;sss*2-1,sss*2-1];%two genes assumed to be distinct: 2 4 6 8 10;1 3 5 7 9

                [xx,n,sigma] = Solver(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax);
                N = 2;
                index = size(n,1);  % The number of the stable states
                alpha = zeros(index,1);  % The weight of the stable states
                mu = zeros(index,N);  % The mean value of the Gaussian density function
                for i = 1:index
                    mu(i,:) = xx(n(i,1),:); 
                    sigma0{i} = reshape(sigma(:,:,i),N,N)';  % The covariance value of the Gaussian density function
                    alpha(i) = n(i,2)/sum(n(:,2));   
                end
                %  mu(find(alpha<=0.01),:)=[]; %disregard states with low weights
                %  sigma0(find(alpha<=0.01))=[];
                %  alpha(find(alpha<=0.01))=[];

                if size(mu,1)==3
                    wei_dist = (alpha-dataweight).^2;    
                    current_loss = 1/2*(dataweight(1)+alpha(1))*1/2*((mu(1,:)-datamu1)*(datasigma(:,:,1)\(mu(1,:)-datamu1)')-log(det(datasigma(:,:,1)\sigma0{1}))+trace(datasigma(:,:,1)\sigma0{1})-N)...
                                  +1/2*(dataweight(2)+alpha(2))*1/2*((mu(2,:)-datamu2)*(datasigma(:,:,2)\(mu(2,:)-datamu2)')-log(det(datasigma(:,:,2)\sigma0{2}))+trace(datasigma(:,:,2)\sigma0{2})-N)...
                                  +1/2*(dataweight(3)+alpha(3))*1/2*((mu(3,:)-datamu3)*(datasigma(:,:,3)\(mu(3,:)-datamu3)')-log(det(datasigma(:,:,3)\sigma0{3}))+trace(datasigma(:,:,3)\sigma0{3})-N)...
                                  +1/2*lambda_penalty*sum(wei_dist);
                              % term 1/2*(dataweight(i)+alpha(i)) can be replaced by dataweight(i)
                    if current_loss < min_loss
                        min_loss = current_loss;
                        best_G = G;
                        best_g = g;
                        best_S = S;
                    end
                end
            end
        end
    end
end

% best_G = [16;16]; best_g = [4;4]; best_S = [4 4;5 5];
