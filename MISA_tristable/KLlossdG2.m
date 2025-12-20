function dlossdG = KLlossdG2(datamu1,datamu2,datamu3,datasigma1,datasigma2,datasigma3,...
    sigma1,sigma2,sigma3,dSigma1dG,dSigma2dG,dSigma3dG,dv1dG,dv2dG,dv3dG,mu11,mu22,mu33,alpha,dataweight)
 
 dlossdG = zeros(size(mu11,2),1);
 dloss1dG = zeros(size(mu11,2),1); dloss2dG = zeros(size(mu11,2),1); dloss3dG = zeros(size(mu11,2),1);

  for i = 1:size(mu11,2)
 det0 = 0;
 det00 = 0;
 det000 = 0;

 dloss1dG(i) = 1/2*(alpha(1)+dataweight(1))*(mu11-datamu1)*(datasigma1\dv1dG(:,i))+1/2*(alpha(2)+dataweight(2))*(mu22-datamu2)*(datasigma2\dv2dG(:,i))+...
     1/2*(alpha(3)+dataweight(3))*(mu33-datamu3)*(datasigma3\dv3dG(:,i));
  for ee = 1:size(mu11,2)
     sigma0 = sigma1;
     sigma0(:,ee) = dSigma1dG(:,ee,i);
     det0 = det0 + det(sigma0);
     sigma00 = sigma2;
     sigma00(:,ee) = dSigma2dG(:,ee,i);
     det00 = det00 + det(sigma00);
     sigma000 = sigma3;
     sigma000(:,ee) = dSigma3dG(:,ee,i);
     det000 = det000 + det(sigma000);
  end
  dloss2dG(i) = -1/2*1/2*(alpha(1)+dataweight(1))*1/det(sigma1)*det0-1/2*1/2*(alpha(2)+dataweight(2))*1/det(sigma2)*det00-...
                 1/2*1/2*(alpha(3)+dataweight(3))*1/det(sigma3)*det000;
  dloss3dG(i) = 1/2*1/2*(alpha(1)+dataweight(1))*sum(sum(inv(datasigma1).*dSigma1dG(:,:,i)'))+...
                1/2*1/2*(alpha(2)+dataweight(2))*sum(sum(inv(datasigma2).*dSigma2dG(:,:,i)'))+...
                1/2*1/2*(alpha(3)+dataweight(3))*sum(sum(inv(datasigma3).*dSigma3dG(:,:,i)'));
  dlossdG(i) = dloss1dG(i)+dloss2dG(i)+dloss3dG(i);
 end
 