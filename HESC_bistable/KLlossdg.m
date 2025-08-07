function dlossdg = KLlossdg(datamu1,datamu2,datasigma1,datasigma2,...
                   sigma1,sigma2,dSigma1dg,dSigma2dg,...
                   dv1dg,dv2dg,mu1,mu2,alpha,dataweight) 

dlossdg = zeros(size(mu1,2),1);
dloss1dg = zeros(size(mu1,2),1); dloss2dg = zeros(size(mu1,2),1); dloss3dg = zeros(size(mu1,2),1);

 for i = 1:size(mu1,2)
 det0 = 0;
 det00 = 0;

 
 dloss1dg(i) = dataweight(1)*(mu1-datamu1)*(datasigma1\dv1dg(:,i))+dataweight(2)*(mu2-datamu2)*(datasigma2\dv2dg(:,i));
 for ee = 1:size(mu1,2)
     sigma0 = sigma1;
     sigma0(:,ee) = dSigma1dg(:,ee,i);
     det0 = det0 + det(sigma0);
     sigma00 = sigma2;
     sigma00(:,ee) = dSigma2dg(:,ee,i);
     det00 = det00 + det(sigma00);
   
  end
  dloss2dg(i) = -1/2*dataweight(1)*1/det(sigma1)*det0-1/2*dataweight(2)*1/det(sigma2)*det00;
  dloss3dg(i) = 1/2*dataweight(1)*sum(sum(inv(datasigma1).*dSigma1dg(:,:,i)'))+...
                1/2*dataweight(2)*sum(sum(inv(datasigma2).*dSigma2dg(:,:,i)'));  
  dlossdg(i) = dloss1dg(i)+dloss2dg(i)+dloss3dg(i);
 end
end
 
 
 
 