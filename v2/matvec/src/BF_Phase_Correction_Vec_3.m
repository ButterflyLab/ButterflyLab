function U = BF_Phase_Correction_Vec_3(U,tau,posst,valst,posst2,valst2,posst3,valst3,md0)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.
if max(abs(mod(U(posst),tau)-mod(valst(:),tau)))>1e-5 && max(abs(mod(U(posst)+tau/2,tau)-mod(valst(:)+tau/2,tau)))>1e-5
    if 0
        mod(U(posst),tau)
        mod(valst,tau)
        error('conflicts');
    else
        U(posst) = valst;
    end
end
if max(abs(mod(U(posst2),tau)-mod(valst2(:),tau)))>1e-5 && max(abs(mod(U(posst2)+tau/2,tau)-mod(valst2(:)+tau/2,tau)))>1e-5
    if 0
        mod(U(posst2),tau)
        mod(valst2,tau)
        error('conflicts2');
    else
        U(posst2) = valst2;
    end
end
if max(abs(mod(U(posst3),tau)-mod(valst3(:),tau)))>1e-5 && max(abs(mod(U(posst3)+tau/2,tau)-mod(valst3(:)+tau/2,tau)))>1e-5
    if 0
        mod(U(posst3),tau)
        mod(valst3,tau)
        error('conflicts3');
    else
        U(posst3) = valst3;
    end
end
U(posst) = valst;
U(posst2) = valst2;
U(posst3) = valst3;
% determing starting values
for cntd = 1:numel(posst)
    mk2nd = md0(cntd);
    der1st = U(posst3(cntd))-U(posst3(cntd)-1);
    num = 0;
    if cntd == numel(posst)
        ed = numel(U);
    else
        ed = posst(cntd+1)-1;
    end
    for cnt = posst3(cntd)+1:ed
        num = num + 1;
        aug = round(((U(cnt-1)+der1st+mk2nd-U(cnt))+tau/2)/tau)*tau;
        vec = aug+(-2:2)*tau+U(cnt);
        [~,pos] = min(abs(vec-U(cnt-1)-der1st-mk2nd));
        U(cnt) = vec(pos);
        % consider second order derivative
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
        der1st = U(cnt)-U(cnt-1);
        mk2nd = der2nd;
    end
end
end


