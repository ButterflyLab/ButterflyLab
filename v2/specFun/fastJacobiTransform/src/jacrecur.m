function [vals] = jacrecur(nts,ts,wghts,n,da,db)

%   Evaluate the functions (1) for dnu=0,1,...,n at a specified point t
%   using the well-known recurrence relations.

nt = zeros(nts,1);
nu = [1:n]';
vals = zeros(nts,n+1);
if  da + db == -1 | da == -1 | db == -1
    vals(:,1) = (cos(ts/2).^(da+1/2)).*(sin(ts/2).^(db+1/2));
    nm = norm(vals(:,1));
    vals(:,1) = vals(:,1)/nm;
    vals(:,1) = vals(:,1)./sqrt(wghts);
else
    C = sqrt((1+da+db)*gamma(1)*gamma(1+da+db)/gamma(1+da)/gamma(1+db));
    vals(:,1) = C*(cos(ts/2).^(da+1/2)).*(sin(ts/2).^(db+1/2));
end
vals(:,2:end) = interpjac1(nt,ts,nu,da,db,-1);
ns = exp(1i*ts*nu');
vals(:,2:end) = real(vals(:,2:end).*ns);

%nu = [0:n]';
%vals = zeros(nts,n+1);
%vals = interpjac1(nt,ts,nu,da,db,-1);
%ns = exp(1i*ts*nu');
%vals = real(vals.*ns);

% if  da == -0.5 & db ==-0.5
%     dnu = n;
%     vals = zeros(nts,n+1);
%     vals(:,1) = ones(nts,1);
%     vals(:,2) = ts;
%     for i = 3:n+1
%         vals(:,i) = 2*ts.*vals(:,i-1)-vals(i-2);
%     end
%     %for i = 1:n+1
%     %    vals(:,i) = vals(:,i)/norm(vals(:,i),2);
%     %end
% else
%     dnu = n;
%     vals = zeros(nts,n+1);
%     for ii =1:nts
%         t    = ts(ii);
%         x    = cos(t);
% 
%         vals(ii,1) = sqrt(2^(-1-da-db)*(1+da+db))*sqrt(gamma(1+da+db)/(gamma(1+da)*gamma(1+db)));
% 
%         if n > 0
%         vals(ii,2) = (sqrt(2^(-1-da-db)*(3+da+db))*(da-db+(2+da+db)*x)*sqrt(gamma(2+da+db)/(gamma(2+da)*gamma(2+db))))/2;
%         end
% 
%         for i=3:n+1
% 
%         dd1 = (sqrt((-1+da+db+2*i)*(1+da+db+2*i))*(4*(-1+da+db)*i*x+4*i^2*x+(da+db)*(da-db+(-2+da+db)*x)))...
%             /(2*sqrt((i*(da+i)*(db+i))/(da+db+i))*(da+db+i)*(-2+da+db+2*i));
% 
%         dd2 =(2^((da+db)/2)*(-1+da+i)*(-1+db+i)*(da+db+2*i))/(i*(da+db+i)*sqrt(-3+da+db+2*i)*(-2+da+db+2*i)...
%             *sqrt((2^(da+db)*(-1+da+i)*(da+i)*(-1+db+i)*(db+i))/((-1+i)*i*(-1+da+db+i)*(da+db+i)*(1+da+db+2*i))));
% 
%         vals(ii,i) = dd1*vals(ii,i-1) - dd2*vals(ii,i-2);
% 
% 
%         end
% 
%         %  Scale by r(t) now
%         rval       = 2^((1+da+db)/2) * cos(t/2)^(db+0.5) * sin(t/2)^(da+0.5);
%         vals(ii,:) = vals(ii,:) * rval;
% 
% 
%     end
%     %for i = 1:n+1
%     %    vals(:,i) = vals(:,i)/norm(vals(:,i),2);
%     %end
% 
end
