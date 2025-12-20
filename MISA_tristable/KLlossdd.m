function dlossdd = KLlossdd(datasigma1,datasigma2,datasigma3,sigma1,sigma2,sigma3,...
    dSigma1dd,dSigma2dd,dSigma3dd,mu1,mu2,mu3,alpha,dataweight)
 dloss1dd = 0;
 det0 = 0;
 det00 = 0;
 det000 = 0;
 
  for ee = 1:size(mu1,2)
     sigma0 = sigma1;
     sigma0(:,ee) = dSigma1dd(:,ee);
     det0 = det0 + det(sigma0);
     sigma00 = sigma2;
     sigma00(:,ee) = dSigma2dd(:,ee);
     det00 = det00 + det(sigma00);
     sigma000 = sigma3;
     sigma000(:,ee) = dSigma3dd(:,ee);
     det000 = det000 + det(sigma000);
  
  end
  dloss2dd = -1/2*1/2*(alpha(1)+dataweight(1))*1/det(sigma1)*det0-1/2*1/2*(alpha(2)+dataweight(2))*1/det(sigma2)*det00...
      -1/2*1/2*(alpha(3)+dataweight(3))*1/det(sigma3)*det000;
  dloss3dd = 1/2*1/2*(alpha(1)+dataweight(1))*sum(sum(inv(datasigma1).*dSigma1dd(:,:)'))+...
             1/2*1/2*(alpha(2)+dataweight(2))*sum(sum(inv(datasigma2).*dSigma2dd(:,:)'))+...
             1/2*1/2*(alpha(3)+dataweight(3))*sum(sum(inv(datasigma3).*dSigma3dd(:,:)'));  
  dlossdd = dloss1dd+dloss2dd+dloss3dd;
 end
 