function Z = HSSBF_Zfun(i,j,num,N,omega,mu0,dl,k,rho,gamma,sc)

% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

i = i(:); j = j(:)';
ind = abs(repmat(i,[1,length(j)])-repmat(j,[length(i),1]));
sz = size(ind);
Z = zeros(sz);

pos = find(ind>=1);
[pos1,pos2] = ind2sub(sz,pos);
%if (abs(i-j)>1) % change >= to >, by Haizhao
Z(pos) = omega*mu0*dl(j(pos2))/4.*besselh(0,2,k*sqrt(sum((rho(:,i(pos1))-rho(:,j(pos2))).^2)))/sc;

pos = find(ind==0);
[pos1,pos2] = ind2sub(sz,pos);
% elseif (abs(i-j)==0)
Z(pos) = omega*mu0*dl(j(pos2))/4.*(1-1i*(2/pi)*log(gamma*k*dl(j(pos2))/4/exp(1.)))/sc;

