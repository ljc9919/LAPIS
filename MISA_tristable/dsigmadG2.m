function [dSigmadG] = dsigmadG2(g,G,k,mu,S,C,D,N,sigma1,dv1dG,Jacob1,dd)
% the last dimension indicates g, the first two indicate sigma
dSigmadG = zeros(N,N,N);
dGdfdx = dG2_dfdx(g,G,mu,S,C,D,N,dv1dG);
for i = 1:N
    A = zeros(N*N);
    b = zeros(N*N,1);
    for pp = 1:N
        oo = pp*N-N+1:pp*N;
        for qq = oo(1):oo(N)
            hh = mod(qq,N);
            if hh == 0
                hh = N;
            end
            A(qq,N*hh-N+1:N*hh) = Jacob1(pp,1:N);
            if hh~=pp
            A(qq,N*pp-N+1:N*pp) = Jacob1(hh,1:N);
            end
        end
    end
     
    for aa = 1:N
    b(aa*N-N+aa) = -sigma1(aa,:)*dGdfdx(aa,:,i)';
        for mm = setdiff(1:N,aa)
            b(aa*N-N+mm) = -sigma1(aa,:)*dGdfdx(mm,:,i)'-sigma1(mm,:)*dGdfdx(aa,:,i)';
        end      
    end

    no = 1; A2 = A;ourrank=zeros(N*(N-1)/2,1);
    for pp = 2:N
        for vv = 1:pp-1
            ourrank(no) = N*(pp-1)+vv; 
            A2(:,N*(vv-1)+pp) = A2(:,ourrank(no))+A2(:,N*(vv-1)+pp);
            A2(:,ourrank(no)) = zeros(size(A2,1),1);no=no+1;
        end
    end
    
    Arank = A2(setdiff(1:size(A,1),ourrank),setdiff(1:size(A,1),ourrank));
    brank = b(setdiff(1:size(A,1),ourrank));
    
    PP = Arank\brank;
    m1=length(PP);
    n2=(-1+sqrt(1+8*m1))/2 ; 
    AA=zeros(n2,n2);
    inde=1;

    for iii=1:n2
        for jjj=iii:n2
        AA(iii,jjj)=PP(inde);
        inde = inde+1;
        end
    end
      AA2 = triu(AA,1);
      AA = AA+AA2';
    
    dSigmadG(:,:,i) = AA;
end
end
