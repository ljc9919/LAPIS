function dlossdS = KLlossds(datamu1,datamu2,datasigma1,datasigma2,...
                   sigma1,sigma2,dSigma1dS,dSigma2dS,...
                   dv1dS,dv2dS,mu1,mu2,alpha,dataweight) 

dlossdS = zeros(size(mu1,2),1);
dloss1dS = zeros(size(mu1,2),1); dloss2dS = zeros(size(mu1,2),1); dloss3dS = zeros(size(mu1,2),1);

 for i = 1:size(mu1,2)
 det0 = 0;
 det00 = 0;

 
 dloss1dS(i) = dataweight(1)*(mu1-datamu1)*(datasigma1\dv1dS(:,i))+dataweight(2)*(mu2-datamu2)*(datasigma2\dv2dS(:,i));

 for ee = 1:size(mu1,2)
     sigma0 = sigma1;
     sigma0(:,ee) = dSigma1dS(:,ee,i);
     det0 = det0 + det(sigma0);
     sigma00 = sigma2;
     sigma00(:,ee) = dSigma2dS(:,ee,i);
     det00 = det00 + det(sigma00);
 
  end
  dloss2dS(i) = -1/2*dataweight(1)*1/det(sigma1)*det0-1/2*dataweight(2)*1/det(sigma2)*det00;
                
  dloss3dS(i) = 1/2*dataweight(1)*sum(sum(inv(datasigma1).*dSigma1dS(:,:,i)'))+...
                1/2*dataweight(2)*sum(sum(inv(datasigma2).*dSigma2dS(:,:,i)'));  
  dlossdS(i) = dloss1dS(i)+dloss2dS(i)+dloss3dS(i);
 end
end
 
 
 
 