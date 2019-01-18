function [Factors,Rcomp] = IBF_multiscale_uniform(fun,xx,kk,NG,tol,disSet,isRecomp)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

% xx and kk are column vectors

if nargin < 6, disSet = []; end;
if nargin < 7, isRecomp = 1; end;

if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end

Dim = size(kk,2);
assert(Dim==1,...
    'xx and kk should be column vectors; only supports 1-dimensional problems in this version.');

Nk = size(kk,1);
Nk = Nk^(1/Dim)*ones(1,Dim);
coronalevels = min(ceil(max(log2(Nk)-6,0)));
kkidglobal = 1:size(kk,1);
kbox(1,:) = min(kk);
kbox(2,:) = max(kk)+(max(kk)-min(kk))./Nk;
kbox = kbox - Nk/2;

Factors = cell(2*coronalevels+1,2);
Rcomp = 0;

for iter = 1:coronalevels
    ckbox = kbox/2^iter;
    
    switch Dim
        case 1
            ck1s = ckbox(1,1)+Nk/2;
            ck1e = ckbox(2,1)+Nk/2;
            
            kkidl = find( kk(:,1)<ck1s );
            kkidr = find( kk(:,1)>=ck1e );
            ckkid = find( kk(:,1)>=ck1s & kk(:,1)<ck1e );
            ssl = 1:numel(kkidl); ssl = ssl(:);
            ssr = 1:numel(kkidr); ssr = ssr(:);
    end
    
    kkcp = kk;
    kkl = kk(kkidl(1),:);
    kkr = kk(kkidr(1),:);
    
    funn = @(t,s) fun(t,s+kkl-1);
    Factors{2*iter-1,2} = kkidglobal(kkidl);
    [Factors{2*iter-1,1},Rcomp_tmp] = IBF_uniform(funn,xx,ssl,NG,tol,disSet,isRecomp,'multiscale');
    Rcomp = Rcomp + Rcomp_tmp*length(kkidl);
    
    funn = @(t,s) fun(t,s+kkr-1);
    Factors{2*iter,2} = kkidglobal(kkidr);
    [Factors{2*iter,1},Rcomp_tmp] = IBF_uniform(funn,xx,ssr,NG,tol,disSet,isRecomp,'multiscale');
    Rcomp = Rcomp + Rcomp_tmp*length(kkidr);
    
    kk = kkcp(ckkid,:);
    kkidglobal = kkidglobal(ckkid);
end

Factors{end,2} = kkidglobal;
Factors{end,1} = fun(xx,kk);
Rcomp = Rcomp + length(kkidglobal);
Rcomp = Rcomp/prod(Nk);

end
