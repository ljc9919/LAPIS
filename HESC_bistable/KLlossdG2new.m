function dlossdG = KLlossdG2new(datamu1,datamu2,datasigma1,datasigma2,...
    sigma1,sigma2,dSigma1dG,dSigma2dG,dv1dG,dv2dG,mu11,mu22,alpha,dataweight,dalphadG,lambda_penalty)
 
 dlossdG = zeros(size(mu11,2),1);
 dloss1dG = zeros(size(mu11,2),1); dloss2dG = zeros(size(mu11,2),1); dloss3dG = zeros(size(mu11,2),1);
 dloss4dG = zeros(size(mu11,2),1); 
 N = size(mu11,2);
 
 for i = 1:size(mu11,2)
 det0 = 0;
 det00 = 0;

 dloss1dG(i) = dataweight(1)*(mu11-datamu1)*(datasigma1\dv1dG(:,i))+dataweight(2)*(mu22-datamu2)*(datasigma2\dv2dG(:,i));
  for ee = 1:size(mu11,2)
     sigma0 = sigma1;
     sigma0(:,ee) = dSigma1dG(:,ee,i);
     det0 = det0 + det(sigma0);
     sigma00 = sigma2;
     sigma00(:,ee) = dSigma2dG(:,ee,i);
     det00 = det00 + det(sigma00);
   
  end
  dloss2dG(i) = -1/2*dataweight(1)*1/det(sigma1)*det0-1/2*dataweight(2)*1/det(sigma2)*det00;
  dloss3dG(i) = 1/2*dataweight(1)*sum(sum(inv(datasigma1).*dSigma1dG(:,:,i)'))+...
                1/2*dataweight(2)*sum(sum(inv(datasigma2).*dSigma2dG(:,:,i)'));
            
  dloss4dG(i) = lambda_penalty*sum((alpha-dataweight).*dalphadG(:,i));
            
  dlossdG(i) = dloss1dG(i)+dloss2dG(i)+dloss3dG(i)+dloss4dG(i);
 end
 