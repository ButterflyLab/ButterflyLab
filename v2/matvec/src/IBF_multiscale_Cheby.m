function [Factors,Rcomp] = IBF_multiscale_Cheby(fun,xx,kk,NG,tol)
% O(N log N) operation and memory complexity.
% The IBF is credited to Yingzhou Li and Haizhao Yang, Interpolative 
% Butterfly Factorization, SIAM J. Sci. Comput., 39(2), A503?A531. 
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang


if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end

Dim = size(kk,2);

Nk = size(kk,1);
Nk = Nk^(1/Dim)*ones(1,Dim);
coronalevels = min(ceil(log2(Nk/16)));
kkidglobal = 1:size(kk,1);
kbox(1,:) = min(kk);
kbox(2,:) = max(kk)+(max(kk)-min(kk))./Nk;

Factors = cell(coronalevels+1,2);
Rcomp = 0;

for iter = 1:coronalevels
    ckbox = kbox/2^iter;
    
    switch Dim
        case 1
            ck1s = ckbox(1,1);
            ck1e = ckbox(2,1);
            
            kkid = find( (kk(:,1)<ck1s | kk(:,1)>=ck1e));
            ckkid = find( kk(:,1)>=ck1s & kk(:,1)<ck1e );
            
        case 2
            ck1s = ckbox(1,1);
            ck1e = ckbox(2,1);
            ck2s = ckbox(1,2);
            ck2e = ckbox(2,2);
            
            kkid = find( (kk(:,1)<ck1s | kk(:,1)>=ck1e) ...
                | (kk(:,2)<ck2s | kk(:,2)>=ck2e) );
            ckkid = find( kk(:,1)>=ck1s & kk(:,1)<ck1e ...
                & kk(:,2)>=ck2s & kk(:,2)<ck2e );
        case 3
            ck1s = ckbox(1,1);
            ck1e = ckbox(2,1);
            ck2s = ckbox(1,2);
            ck2e = ckbox(2,2);
            ck3s = ckbox(1,3);
            ck3e = ckbox(2,3);
            
            kkid = find( (kk(:,1)<ck1s | kk(:,1)>=ck1e) ...
                | (kk(:,2)<ck2s | kk(:,2)>=ck2e) ...
                | (kk(:,3)<ck3s | kk(:,3)>=ck3e) );
            ckkid = find( kk(:,1)>=ck1s & kk(:,1)<ck1e ...
                & kk(:,2)>=ck2s & kk(:,2)<ck2e ...
                & kk(:,3)>=ck3s & kk(:,3)<ck3e);
    end
    
    kkcp = kk;
    kk = kk(kkid,:);
    
    Factors{iter,2} = kkidglobal(kkid);
    [Factors{iter,1},Rcomp_tmp] = IBF_Cheby(fun,xx,kk,NG,tol,'multiscale');
    Rcomp = Rcomp + Rcomp_tmp*length(kkid);
    
    kk = kkcp(ckkid,:);
    kkidglobal = kkidglobal(ckkid);
end

Factors{end,2} = kkidglobal;
Factors{end,1} = fun(xx,kk);
Rcomp = Rcomp + length(kkidglobal);
Rcomp = Rcomp/prod(Nk);

end
