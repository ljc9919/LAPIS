function dlossdd = KLlossdd(datasigma1,datasigma2,sigma1,sigma2,...
    dSigma1dd,dSigma2dd,mu1,mu2,alpha,dataweight)
 dloss1dd = 0;
 det0 = 0;
 det00 = 0;
 
  for ee = 1:size(mu1,2)
     sigma0 = sigma1;
     sigma0(:,ee) = dSigma1dd(:,ee);
     det0 = det0 + det(sigma0);
     sigma00 = sigma2;
     sigma00(:,ee) = dSigma2dd(:,ee);
     det00 = det00 + det(sigma00);
    
  end
  dloss2dd = -1/2*dataweight(1)*1/det(sigma1)*det0-1/2*dataweight(2)*1/det(sigma2)*det00;
      
  dloss3dd = 1/2*dataweight(1)*sum(sum(inv(datasigma1).*dSigma1dd(:,:)'))+...
             1/2*dataweight(2)*sum(sum(inv(datasigma2).*dSigma2dd(:,:)'));  
  dlossdd = dloss1dd+dloss2dd+dloss3dd;
end
 