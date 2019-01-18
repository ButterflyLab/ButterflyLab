function [Factors,Rcomp] = fastMBF(fun,xx,kk,NG,tol)

Dim = size(kk,2);
assert(Dim==2||Dim==3,...
    'fastMBF only supports 2- and 3-dimensional problems.');

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
    
    if Dim == 2
        ck1s = ckbox(1,1);
        ck1e = ckbox(2,1);
        ck2s = ckbox(1,2);
        ck2e = ckbox(2,2);
        
        kkid = find( (kk(:,1)<ck1s | kk(:,1)>=ck1e) ...
            | (kk(:,2)<ck2s | kk(:,2)>=ck2e) );
        ckkid = find( kk(:,1)>=ck1s & kk(:,1)<ck1e ...
            & kk(:,2)>=ck2s & kk(:,2)<ck2e );
    end
    
    if Dim == 3
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
    [Factors{iter,1},Rcomp_tmp] = fastBF(fun,xx,kk,NG,tol,'multiscale');
    Rcomp = Rcomp + Rcomp_tmp*length(kkid);
    
    kk = kkcp(ckkid,:);
    kkidglobal = kkidglobal(ckkid);
end

Factors{end,2} = kkidglobal;
Factors{end,1} = fun(xx,kk);
Rcomp = Rcomp + length(kkidglobal);
Rcomp = Rcomp/prod(Nk);

end
