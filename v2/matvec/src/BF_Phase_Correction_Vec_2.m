function U = BF_Phase_Correction_Vec_2(U,tau,posst,valst)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if max(abs(mod(U(posst),tau)-mod(valst,tau)))>1e-5 & max(abs(mod(U(posst)+tau/2,tau)-mod(valst+tau/2,tau)))>1e-5
    if 0
        mod(U(posst),tau)
        mod(valst,tau)
        error('conflicts');
    else
        U(posst) = valst;
    end
end
U(posst) = valst;
for cntd = 1:numel(posst)
    mk2nd = 0;
    num = 0;
    if abs(cntd-numel(posst))<1e-10
        edPoint = numel(U);
    else
        edPoint = posst(cntd+1)-1;
    end
    der1st = 0;
    for cnt = posst(cntd)+1:edPoint
        num = num + 1;
        aug = round(((U(cnt-1)+der1st+mk2nd-U(cnt))+tau/2)/tau)*tau;
        vec = aug+(-2:2)*tau+U(cnt);
        if num > 1
            if num > 2
                [~,pos] = min(abs(vec-U(cnt-1)-der1st-mk2nd));
                U(cnt) = vec(pos);
                der2nd = (U(cnt)-2*U(cnt-1)+U(cnt-2));
                while abs(der2nd-mk2nd)>0.95*tau
                    vec(pos) = [];
                    if numel(vec) < 1e-5
                        error('empty');
                    end
                    [~,pos] = min(abs(vec-U(cnt-1)-der1st-mk2nd));
                    U(cnt) = vec(pos);
                    der2nd = (U(cnt)-2*U(cnt-1)+U(cnt-2));
                end
            else
                [~,pos] = min(abs((vec-U(cnt-1))-(U(cnt-1)-U(cnt-2))));
                U(cnt) = vec(pos);
                der2nd = (U(cnt)-2*U(cnt-1)+U(cnt-2));
            end
            mk2nd = der2nd;
        else
            [~,pos] = min(abs(vec-U(cnt-1)-der1st-mk2nd));
            U(cnt) = vec(pos);
        end
        der1st = U(cnt)-U(cnt-1);
    end
end