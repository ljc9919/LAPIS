function [xx,n,sigma]=Solver(cycle_index,k,G,g,S,IM,Hilln,AB,C,D,dd,wemax)
N = 2;  %%the dimension of the system
xx = zeros(cycle_index,N);
rng('default')
%%Solve odes from different initial values
for i = 1:cycle_index
    x0 = rand(1,N).*wemax;
    [t,x] = ode15s(@(t,x)force(t,x,k,G,g,S,IM,Hilln,AB),[0,100],x0);
    newx = x(end,:);
    x = inf*ones(1,N);
    while norm( x(end,:)-newx(end,:) ,2 )>1e-7
        x = newx;
        [t,newx] = ode15s(@(t,x)force(t,x,k,G,g,S,IM,Hilln,AB),[0,100],x(end,:));
    end
    xx(i,:) = newx(end,:);
end

%%Finding the stable points
 for q = 1:(cycle_index-1)
     for p = (q+1):cycle_index
         if norm(xx(q,:)-xx(p,:),'fro')<10^-1
             xx(p,:) = xx(q,:);
         end
     end
 end
stable_point = unique(xx(:,:),'rows');

sigma = [];
for i = 1:size(stable_point,1)
    [m] = find(xx(:,2)==stable_point(i,2));
    if length(m)>=1
        disp(strcat(num2str(stable_point(i,:)),' repeat ',num2str(length(m)),' times',' the location in the row xx is' ,mat2str(m)))
    end
   n(i,1) = m(1);
   n(i,2) = length(m);
   sigma(:,:,i) = Sigma(xx(m(1),:)',IM,C,D,S,Hilln,G,g,k,dd);
end

%Arrange the index n from large to small by the first element
m = size(n,1);
if(m ~= 1)
    for i = 1:m
        tran = n(i,:);
        flag = i;
        if( i ~= m )
            for j = i+1:m
                if( xx(n(j,1),1) > xx(tran(1),1) )
                    tran = n(j,:);
                    flag = j;
                end
            end
        end
        n(flag,:) = n(i,:);
        n(i,:) = tran;
    end
end
%%stable state
SS = zeros(size(xx,2),size(n,1));
for i = 1:size(n,1)
    SS(:,i) = xx(n(i,1),:)';
end
end