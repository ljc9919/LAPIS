f size(mu,1)~=2
       continue;
    else
        alphaold = alpha; Gold = G; gold = g; Sold = S; ddold = dd;
        [loss,dlossdg,dlossdG,dlossdS,dlossdd] = Solveloss0(dataweight,datamu1,datamu2,datasigma,alpha,mu,sigma,g,G,k,S,C,D,IM,N,dd,Hilln,lambda_penalty);
        
        Q = 1; Cval = loss; tau1 = 1e-4; tau2 = 1e-4; tau3 = 1e-4; tau4 = 1e-5;
        for iter = 1 : itermax
            xp = [g;G;S(1:5,1);dd];     fp = loss;     gp = [dlossdg;dlossdG;dlossdS(:);dlossdd];
            nls = 1;
            gold = xp(1:N); Gold = xp(N+1:2*N); Sold = repmat(xp(2*N+1:3*N),1,5); ddold = xp(end);
            Cval = loss;
            while 1
                x = xp-[tau1*ones(N,1);tau2*ones(N,1);tau3*ones(N,1);tau4].*gp;
                g = x(1:N); G = x(N+1:2*N); S = repmat(x(2*N+1:3*N),1,5); dd = x(end);
                [mu,sigma,alpha] = Solveode(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax);
                if nls>=10
                    break
                end
                if size(mu,1)==2
                    [loss,dlossdg,dlossdG,dlossdS,dlossdd] = Solveloss(dataweight,datamu1,datamu2,datasigma,alpha,alphaold,mu,sigma,g,gold,G,Gold,k,S,Sold,C,D,IM,N,dd,ddold,Hilln,lambda_penalty);
                    nrmG1  = norm(dlossdg); nrmG2 = norm(dlossdG); nrmG3  = norm(dlossdS); nrmG4  = norm(dlossdd);
                    if loss <= Cval - (tau1*rhols*nrmG1^2+tau2*rhols*nrmG2^2+tau3*rhols*nrmG3^2+tau4*rhols*nrmG4^2) || nls >= 10
                        break
                    end
                    tau = eta*[tau1,tau2,tau3,tau4];
                    tau1 = tau(1); tau2 = tau(2); tau3 = tau(3); tau4 = tau(4);
                    nls = nls+1;
                else
                    tau = eta*[tau1,tau2,tau3,tau4];
                    tau1 = tau(1); tau2 = tau(2); tau3 = tau(3); tau4 = tau(4);
                    nls = nls+1;
                end
            end

            if size(mu,1)~=2
                break;
            end
    
          	nrmG  = norm([dlossdg;dlossdG;dlossdS;dlossdd], 'fro');
            s = x - xp; XDiff = norm(s);
            FDiff = abs(fp-loss)/(abs(fp)+1);
    
            if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol
                out.msg = 'converge';
                break;
            end
           	y = [dlossdg;dlossdG;dlossdS;dlossdd] - gp;
            sy1 = abs(dot(s(1:N),y(1:N)));    %tau = opts.tau;
            sy2 = abs(dot(s(N+1:2*N),y(N+1:2*N)));
            sy3 = abs(dot(s(2*N+1:3*N),y(2*N+1:3*N)));
            sy4 = abs(dot(s(3*N+1),y(3*N+1)));
                
           if sy1 > 0
                if mod(iter,2)==0; tau1 = abs(sum(sum(s(1:N).*s(1:N))))/sy1;
                else
                    tau1 = sy1/abs(sum(sum(y(1:N).*y(1:N)))); 
                end
           end
     
            if sy2 > 0
                 if mod(iter,2)==0; tau2 = abs(sum(sum(s(N+1:2*N).*s(N+1:2*N))))/sy2;
                 else
                     tau2 = sy2/abs(sum(sum(y(N+1:2*N).*y(N+1:2*N)))); 
                 end
            end

            if sy3 > 0
                  if mod(iter,2)==0; tau3 = abs(sum(sum(s(2*N+1:3*N).*s(2*N+1:3*N))))/sy3;
                  else
                      tau3 = sy3/abs(sum(sum(y(2*N+1:3*N).*y(2*N+1:3*N))));
                  end
            end

            if sy4 > 0
                  if mod(iter,2)==0; tau4 = abs(sum(sum(s(3*N+1).*s(3*N+1))))/sy4;
                  else
                      tau4 = sy4/abs(sum(sum(y(3*N+1).*y(3*N+1))));
                  end
            end


            tau1 = max(min(tau1, 1e-4), 1e-7);
            tau2 = max(min(tau2, 1e-4), 1e-7);
            tau3 = max(min(tau3, 1e-4), 1e-7);
            tau4 = max(min(tau4, 1e-5), 1e-7);

            Qp = Q; Q = gamma*Qp + 1; 
            lambda_penalty = shrinkage*lambda_penalty;
            lovepara{iter,sample} = x; loveloss{iter,sample} = loss; loveweight{iter,sample} = alpha;
            lovenormG{iter,sample} = nrmG; lovestep{iter,sample} = [tau1,tau2,tau3,tau4]; lovemu{iter,sample} = mu;
            if (iter >= 1)
                fprintf('%4d \t %7.6e \t %3.2e \t %3.2e \t %3.2e \t %2d\n', ...
                            iter, loss, nrmG, XDiff, FDiff, nls);
            end

            if Cval<loss-5 || loss<0
                break;
            end
        end
	end
 end

save second_refining_ND.mat