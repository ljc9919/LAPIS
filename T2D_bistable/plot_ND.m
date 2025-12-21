cycle_index = 1000;
g = [0.1989 
0.4183 
0.3172 
0.2927 
0.0934 
];
G = [0.6523 
0.8605 
0.6549 
0.5691 
0.5468 
];
k = ones(5,1);
S = repmat([0.2931 
0.2728 
0.2962 
0.2400 
0.2549 
],1,5);

dd = 0.0265;

IM=sign([0	0.604328586	0	0	0	0.218094675	0.335669628	0	0	0	0
0.604328586	0	-0.448255943	0	0	0.320282934	0.505490034	0	-0.401133123	0	0
0	-0.448255943	0	0.16837281	0.553556606	0	0	-0.532746333	0	0.748115209	0
0	0	0.16837281	0	0.335283464	0	0.337270073	0	-0.509976426	0	0
0	0	0.553556606	0.335283464	0	0	0	-0.581706519	0	0.805509325	0
0.218094675	0.320282934	0	0	0	0	0.492959914	0.596376723	-0.45860916	0	0
0	0	0	0.337270073	0	0.492959914	0	0	-0.488976146	0	0
0	0	-0.532746333	0	-0.581706519	0.596376723	0	0	0	-0.5090657	0
0	0	0	-0.509976426	0	-0.45860916	0	0	0	0	0
0	0	0.748115209	0	0.805509325	0	0	-0.5090657	0	0	0
0	0	0	0	0	0	0	0	0.567017303	0	0]);
IM=IM(1:5,1:5);

Hilln = 4*ones(N,N);

AB1 = [0	0.38329688	0.186955249	0	0	0	0	0	0	0	0	0.049311856
0.38329688	0	0.604328586	0	0	0	0	0.093643101	0	0	0	0
0.186955249	0.604328586	0	-0.448255943	0	0	0.221507797	0.120112039	0	0	0	0
0	0	-0.448255943	0	0.16837281	0.553556606	0	0	0	0	0.370807464	0
0	0	0	0.16837281	0	0.335283464	0	0	0	0	0	0
0	0	0	0.553556606	0.335283464	0	0	0	-0.178603499	0	0.425933659	0.175174329
0	0	0.221507797	0	0	0	0	0	0.154761295	-0.032780494	0	0
0	0	0.120112039	0	0	0	0	0	0	-0.227451508	0	0
0	0	0	0	0	-0.178603499	0.154761295	0	0	-0.010797602	-0.171342499	0
0	0	0	0	0	0	-0.032780494	-0.227451508	-0.010797602	0	0.073152125	0.116332668
0	0	0	0.370807464	0	0.425933659	0	0	-0.171342499	0.073152125	0	0.188316026
0.049311856	0	0	0	0	0.175174329	0	0	0	0.116332668	0.188316026	0]; 
AB = AB1; AB(find(AB1<0)) = -AB(find(AB1<0));
AB = AB(2:6,2:6);AB1 = AB1(2:6,2:6);

    cno = (AB1(:)>0);
    C = zeros(N,N);
    C(cno) = AB1(cno);

    dno = (AB1(:)<0);
    D = zeros(N,N);
    D(dno) = -AB1(dno);
    
    wemax  = 1;
[xx,n,sigma,ycell,action,ActionVal]=Solver5(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax);

for i=1:index
%The mean value of each stable state
    mu(i,:)=xx(n(i,1),:); 
    sigma0{i}=reshape(sigma(:,:,i),N,N)';
%The weight of each stable state
    alpha(i)=n(i,2)/sum(n(:,2));   %%%%%°??a??alpha?ó??è￥￡?muoíalphaê???ó|μ??￡
end

   V = Vtrue;
   sigma0_pca=cell(index,1);
    mu_pca=zeros(index,2);
    for i=1:index
        mu_pca(i,:)=V'*mu(i,:)';
        sigma0_pca{i}=V'*sigma0{i}*V;
    end
    
    %% plot the landscape
    %y_max=[2,1.5]; %% Range of the landscape
    %y_min=[-1,-1.5];
    y_max=[2,2]; %% Range of the landscape
    y_min=[-1,-2];
    step=(y_max-y_min)/200; %% Length of the step
    [a1,a2]=meshgrid(y_min(1):step(1):y_max(1),y_min(2):step(2):y_max(2)); %% Grid
    [s1,s2]=size(a1);
    P=zeros(s1,s2);
    z=zeros(s1,s2);
    for kk=1:index
        sig=sigma0_pca{kk};
        x_wen=mu_pca(kk,:);
        for i=1:s1
            for j=1:s2
                z(i,j)=multivariate_normal_distribution([a1(i,j);a2(i,j)],x_wen',sig,2);  %% Normal distribution
            end
        end

        P=P+z*alpha(kk);
    end
    P=real(P);
    P=P/sum(sum(P));
 
    surf(a1,a2,-log(max(P,10^-8)));   %% Plot landscape
    shading interp
%     xlabel('HNF1A')
%     ylabel('TCF4')
    xlabel('PC1')
    ylabel('PC2')
    zlabel('U')
    %axis([-1 1.5 -1 1 6 17])
    
        y12=V'*ycell{1,2};
        y21=V'*ycell{2,1};
        view([-25,75])
        hold on
        z3path=griddata(a1,a2,-log(max(P,10^-8)),y12(1,:),y12(2,:));
        plot3(y12(1,:),y12(2,:),z3path+0.5,'w','LineWidth',2);

        z3path=griddata(a1,a2,-log(max(P,10^-8)),y21(1,:),y21(2,:));
        plot3(y21(1,:),y21(2,:),z3path+0.5,'Color',[0.85,0.43,0.83],'LineWidth',2);

        view([-25 75])
        set(gcf,'outerposition', [100 100 800 650]);

        
        %%
           plotycell=ycell{1,2};
    for i=1:5
        yy(i,:)=(plotycell(i,:)-min(plotycell(i,:)))./(max(plotycell(i,:))-min(plotycell(i,:)));
    end
    
    heatmap(yy)
    
    for i=1:5
        for j=1:30
    zz(i,2*j-1)=yy(i,j);zz(i,2*j)=(yy(i,j)+yy(i,j+1))/2;
        end
    end
    
    h=heatmap(zz(:,1:60))
    h.YDisplayLabels = {'HNF4A', 'HNF1A', 'TCF4', 'NEUROD1', 'NFIA'};
    h.XDisplayLabels = repmat({''}, 1, 60);
