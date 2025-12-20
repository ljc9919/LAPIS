clear,clc
tic()
mutrue1 = [16.9089	0.0151;11.7317	12.7648;0.0637	14.8343]; %(highlow,highhigh,lowhigh)ground_truth mu value
mutrue = reshape(mutrue1,1,6);
datamu1 = mutrue1(1,:);   datamu2 = mutrue1(2,:);   datamu3 = mutrue1(3,:); 

dataweight = [0.35500 0.38567 0.25933]';%ground_truth weight value
datasigma(:,:,1) = [0.4000	-0.0036;-0.0036	0.4184];
datasigma(:,:,2) = [0.4585	-0.0141;-0.0141	0.4406];
datasigma(:,:,3) = [0.4087	-0.0007;-0.0007	0.4000];%ground_truth sigma value

N = 2;%Dim
IM = [1 -1;-1 1]; %prior network structure
Hilln = 4*ones(N,N); %Hill coefficient
AB = IM; AB(find(IM<0)) = -AB(find(IM<0));
cno = (IM(:)>0); C = zeros(N,N); C(cno) = IM(cno);%activation structure
dno = (IM(:)<0); D = zeros(N,N); D(dno) = -IM(dno);%inhibition structure

itermax = 1000;
k = 1*ones(N,1);
rng(1);
randno = randn(200,1);
cycle_index = 3000;

for sample=15:20  %1:20 results in longer computation time
    G = [16+randno(sample)*2/3 16+randno(sample+50)*2/3]'; % Generated random initial parameters for each trial
    g = [4+randno(sample+100)*1/3 4+randno(sample+150)*1/3]';
    g(find(g<=0)) = 0; % non-negativity constraint
    g(find(g>G)) = G(find(g>G)); % upper bound constraint
    S = [4 4;5 5]; % can be perturbed, here fixed to avoid phase transitions

    wemax = 10*ones(1,N);
    lambda_penalty = 0; % hyperparameters, more accurate fit to parameters, rather than weights

    no = 1; % iter
    dd = 0.1; % can also be perturbed

    for iter = 1:itermax
        [xx,n,sigma]=Solver(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax);

        N = 2;
        index=size(n,1);  %% The number of the stable states
        alpha=zeros(index,1);  %% The weight of the stable states
        mu=zeros(index,N);  %% The mean value of the Gaussian density function
        
        for i=1:index
            mu(i,:)=xx(n(i,1),:); 
            sigma0{i}=reshape(sigma(:,:,i),N,N)';
            alpha(i)=n(i,2)/sum(n(:,2));   
        end
        %      mu(find(alpha<=0.01),:)=[];
        %      sigma0(find(alpha<=0.01))=[];
        %      alpha(find(alpha<=0.01))=[];
        
        if size(mu,1)==3 %Only trials where the state count matched the target
            wei_dist = (alpha-dataweight).^2;

            loss = 1/2*(dataweight(1)+alpha(1))*1/2*((mu(1,:)-datamu1)*(datasigma(:,:,1)\(mu(1,:)-datamu1)')-log(det(datasigma(:,:,1)\sigma0{1}))+trace(datasigma(:,:,1)\sigma0{1})-N)...
                  +1/2*(dataweight(2)+alpha(2))*1/2*((mu(2,:)-datamu2)*(datasigma(:,:,2)\(mu(2,:)-datamu2)')-log(det(datasigma(:,:,2)\sigma0{2}))+trace(datasigma(:,:,2)\sigma0{2})-N)...
                  +1/2*(dataweight(3)+alpha(3))*1/2*((mu(3,:)-datamu3)*(datasigma(:,:,3)\(mu(3,:)-datamu3)')-log(det(datasigma(:,:,3)\sigma0{3}))+trace(datasigma(:,:,3)\sigma0{3})-N)...
                  +1/2*lambda_penalty*sum(wei_dist); % loss function
                  
            [dv1dg,dv2dg,dv3dg] = dvdg(g,G,k,mu(1,:),mu(2,:),mu(3,:),S,C,D,N); % explicitly compute the gradient of the mean with respect to g, G, and S
            [dv1dG,dv2dG,dv3dG] = dvdG2(g,G,k,mu(1,:),mu(2,:),mu(3,:),S,C,D,N);
            [dv1ds,dv2ds,dv3ds] = dvds(g,G,k,mu(1,:),mu(2,:),mu(3,:),S,C,D,N);

            [Jacob1]=Jacobi(mu(1,:),IM,C,D,S,Hilln,G,g,k);
            [Jacob2]=Jacobi(mu(2,:),IM,C,D,S,Hilln,G,g,k);
            [Jacob3]=Jacobi(mu(3,:),IM,C,D,S,Hilln,G,g,k);

            dSigma1dg=dsigmadg(g,G,k,mu(1,:),S,C,D,N,sigma(:,:,1),dv1dg,Jacob1,dd); % explicitly compute the gradient of the covariance with respect to g, G, S and d
            dSigma2dg=dsigmadg(g,G,k,mu(2,:),S,C,D,N,sigma(:,:,2),dv2dg,Jacob2,dd);
            dSigma3dg=dsigmadg(g,G,k,mu(3,:),S,C,D,N,sigma(:,:,3),dv3dg,Jacob3,dd); 

            dSigma1dG=dsigmadG2(g,G,k,mu(1,:),S,C,D,N,sigma(:,:,1),dv1dG,Jacob1,dd); 
            dSigma2dG=dsigmadG2(g,G,k,mu(2,:),S,C,D,N,sigma(:,:,2),dv2dG,Jacob2,dd);
            dSigma3dG=dsigmadG2(g,G,k,mu(3,:),S,C,D,N,sigma(:,:,3),dv3dG,Jacob3,dd); 

            dSigma1dS=dsigmads(g,G,k,mu(1,:),S,C,D,N,sigma(:,:,1),dv1ds,Jacob1,dd);  
            dSigma2dS=dsigmads(g,G,k,mu(2,:),S,C,D,N,sigma(:,:,2),dv2ds,Jacob2,dd);
            dSigma3dS=dsigmads(g,G,k,mu(3,:),S,C,D,N,sigma(:,:,3),dv3ds,Jacob3,dd);  

            dSigma1dd=dsigmadd(N,Jacob1); dSigma2dd=dsigmadd(N,Jacob2);
            dSigma3dd=dsigmadd(N,Jacob3); 

            dlossdg = KLlossdg(datamu1,datamu2,datamu3,datasigma(:,:,1),datasigma(:,:,2),datasigma(:,:,3),...  % explicitly compute the gradient of the loss with respect to g, G, S and d
                sigma(:,:,1),sigma(:,:,2),sigma(:,:,3),dSigma1dg,dSigma2dg,dSigma3dg,...
                dv1dg,dv2dg,dv3dg,mu(1,:),mu(2,:),mu(3,:),alpha,dataweight); 

            dlossdG = KLlossdG2(datamu1,datamu2,datamu3,datasigma(:,:,1),datasigma(:,:,2),datasigma(:,:,3),...
                sigma(:,:,1),sigma(:,:,2),sigma(:,:,3),dSigma1dG,dSigma2dG,dSigma3dG,...
                dv1dG,dv2dG,dv3dG,mu(1,:),mu(2,:),mu(3,:),alpha,dataweight); 

            dlossdS = KLlossds(datamu1,datamu2,datamu3,datasigma(:,:,1),datasigma(:,:,2),datasigma(:,:,3),...
                sigma(:,:,1),sigma(:,:,2),sigma(:,:,3),dSigma1dS,dSigma2dS,dSigma3dS,...
                dv1ds,dv2ds,dv3ds,mu(1,:),mu(2,:),mu(3,:),alpha,dataweight); 

            dlossdd = KLlossdd(datasigma(:,:,1),datasigma(:,:,2),datasigma(:,:,3),sigma(:,:,1),sigma(:,:,2),sigma(:,:,3),...
                dSigma1dd,dSigma2dd,dSigma3dd,mu(1,:),mu(2,:),mu(3,:),alpha,dataweight);

            if iter<20
            lambdag = 10^(floor(log10(max(abs(dlossdg(:)))))-2.5); % adaptive step size
            lambdaG = 10^(floor(log10(max(abs(dlossdG(:)))))-2.5); 
            lambdaS = 10^(floor(log10(max(abs(dlossdS(:)))))-2.5); 
            lambdad = 10^(floor(log10(max(abs(dlossdd(:)))))-2.5);

            else
            lambdag = 10^(-floor(log10(max(abs(dlossdg(:)))))-2.5); % adaptive step size
            lambdaG = 10^(-floor(log10(max(abs(dlossdG(:)))))-2.5); 
            lambdaS = 10^(-floor(log10(max(abs(dlossdS(:)))))-2.5); 
            lambdad = 10^(-floor(log10(max(abs(dlossdd(:)))))-2.5);
            end

             g1 = max(g - lambdag*dlossdg,zeros(N,1));  % gradient descent 
             G1 = max(G - lambdaG*dlossdG,zeros(N,1));
             S1_vec = max(S(:,1) - lambdaS*dlossdS,zeros(N,1)); S1 = S1_vec*ones(1,N);
             dd1 = max(dd - lambdad*dlossdd,0);

             g = g1; G = G1; dd = dd1; S = S1;
             ourloss(no,sample) = loss; % gradient descent 
             ourmu{no,sample} = mu; oursigma{no,sample} = sigma; 
             ourg{no,sample} = g; ourG{no,sample} = G; ourd(no,sample) = dd; ourS{no,sample} = S; % recorded the results of every iteration across all trials
             no = no+1;
        else
             sample = sample+1;break;
        end
    end
end
