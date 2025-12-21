function [loss,dlossdg,dlossdG,dlossdS,dlossdd]=Solveloss0(dataweight,datamu1,datamu2,datasigma,alpha,mu,sigma,g,G,k,S,C,D,IM,N,dd,Hilln,lambda_penalty)
    wei_dist = (alpha-dataweight).^2;
    
    loss = dataweight(1)*1/2*((mu(1,:)-datamu1)*(datasigma(:,:,1)\(mu(1,:)-datamu1)')-log(det(datasigma(:,:,1)\sigma(:,:,1)))+trace(datasigma(:,:,1)\sigma(:,:,1))-N)...
          +dataweight(2)*1/2*((mu(2,:)-datamu2)*(datasigma(:,:,2)\(mu(2,:)-datamu2)')-log(det(datasigma(:,:,2)\sigma(:,:,2)))+trace(datasigma(:,:,2)\sigma(:,:,2))-N)...
          +1/2*lambda_penalty*sum(wei_dist);
    
    [dv1dg,dv2dg] = dvdg(g,G,k,mu(1,:),mu(2,:),S,C/max(sum(C)),D/max(sum(D)),N);  
    [dv1dG,dv2dG] = dvdG2(g,G,k,mu(1,:),mu(2,:),S,C/max(sum(C)),D/max(sum(D)),N);
    [dv1ds,dv2ds] = dvds(g,G,k,mu(1,:),mu(2,:),S,C/max(sum(C)),D/max(sum(D)),N);

    [Jacob1]=Jacobi(mu(1,:),IM,C/max(sum(C)),D/max(sum(D)),S,Hilln,G,g,k);
    [Jacob2]=Jacobi(mu(2,:),IM,C/max(sum(C)),D/max(sum(D)),S,Hilln,G,g,k);
    
    dSigma1dg=dsigmadg(g,G,k,mu(1,:),S,C/max(sum(C)),D/max(sum(D)),N,sigma(:,:,1),dv1dg,Jacob1,dd); dSigma2dg=dsigmadg(g,G,k,mu(2,:),S,C/max(sum(C)),D/max(sum(D)),N,sigma(:,:,2),dv2dg,Jacob2,dd);
   
    dSigma1dG=dsigmadG2(g,G,k,mu(1,:),S,C/max(sum(C)),D/max(sum(D)),N,sigma(:,:,1),dv1dG,Jacob1,dd); dSigma2dG=dsigmadG2(g,G,k,mu(2,:),S,C/max(sum(C)),D/max(sum(D)),N,sigma(:,:,2),dv2dG,Jacob2,dd);
    
    dSigma1dS=dsigmads(g,G,k,mu(1,:),S,C/max(sum(C)),D/max(sum(D)),N,sigma(:,:,1),dv1ds,Jacob1,dd);  dSigma2dS=dsigmads(g,G,k,mu(2,:),S,C/max(sum(C)),D/max(sum(D)),N,sigma(:,:,2),dv2ds,Jacob2,dd);
     
    dSigma1dd=dsigmadd(N,Jacob1); dSigma2dd=dsigmadd(N,Jacob2);
    
    dlossdg = KLlossdg(datamu1,datamu2,datasigma(:,:,1),datasigma(:,:,2),...
        sigma(:,:,1),sigma(:,:,2),dSigma1dg,dSigma2dg,...
        dv1dg,dv2dg,mu(1,:),mu(2,:),alpha,dataweight); 
    
    dlossdG = KLlossdG2(datamu1,datamu2,datasigma(:,:,1),datasigma(:,:,2),...
        sigma(:,:,1),sigma(:,:,2),dSigma1dG,dSigma2dG,...
        dv1dG,dv2dG,mu(1,:),mu(2,:),alpha,dataweight); 
    
    dlossdS = KLlossds(datamu1,datamu2,datasigma(:,:,1),datasigma(:,:,2),...
        sigma(:,:,1),sigma(:,:,2),dSigma1dS,dSigma2dS,...
        dv1ds,dv2ds,mu(1,:),mu(2,:),alpha,dataweight); 
 
    dlossdd = KLlossdd(datasigma(:,:,1),datasigma(:,:,2),sigma(:,:,1),sigma(:,:,2),...
        dSigma1dd,dSigma2dd,mu(1,:),mu(2,:),alpha,dataweight);
end
    