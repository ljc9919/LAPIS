function dlossdS = KLlossds(datamu1,datamu2,datamu3,datasigma1,datasigma2,datasigma3,...
                   sigma1,sigma2,sigma3,dSigma1dS,dSigma2dS,dSigma3dS,...
                   dv1dS,dv2dS,dv3dS,mu1,mu2,mu3,alpha,dataweight) 

dlossdS = zeros(size(mu1,2),1);
dloss1dS = zeros(size(mu1,2),1); dloss2dS = zeros(size(mu1,2),1); dloss3dS = zeros(size(mu1,2),1);

 for i = 1:size(mu1,2)
 det0 = 0;
 det00 = 0;
 det000 = 0;
 
 dloss1dS(i) = 1/2*(alpha(1)+dataweight(1))*(mu1-datamu1)*(datasigma1\dv1dS(:,i))+1/2*(alpha(2)+dataweight(2))*(mu2-datamu2)*(datasigma2\dv2dS(:,i))...
     +1/2*(alpha(3)+dataweight(3))*(mu3-datamu3)*(datasigma3\dv3dS(:,i));

 for ee = 1:size(mu1,2)
     sigma0 = sigma1;
     sigma0(:,ee) = dSigma1dS(:,ee,i);
     det0 = det0 + det(sigma0);
     sigma00 = sigma2;
     sigma00(:,ee) = dSigma2dS(:,ee,i);
     det00 = det00 + det(sigma00);
     sigma000 = sigma3;
     sigma000(:,ee) = dSigma3dS(:,ee,i);
     det000 = det000 + det(sigma000);
    
  end
  dloss2dS(i) = -1/2*1/2*(alpha(1)+dataweight(1))*1/det(sigma1)*det0-1/2*1/2*(alpha(2)+dataweight(2))*1/det(sigma2)*det00...
                -1/2*1/2*(alpha(3)+dataweight(3))*1/det(sigma3)*det000;
  dloss3dS(i) = 1/2*1/2*(alpha(1)+dataweight(1))*sum(sum(inv(datasigma1).*dSigma1dS(:,:,i)'))+...
                1/2*1/2*(alpha(2)+dataweight(2))*sum(sum(inv(datasigma2).*dSigma2dS(:,:,i)'))+...
                1/2*1/2*(alpha(3)+dataweight(3))*sum(sum(inv(datasigma3).*dSigma3dS(:,:,i)'));  
  dlossdS(i) = dloss1dS(i)+dloss2dS(i)+dloss3dS(i);
 end
end
 
 
 
 