function [U,disPosu] = BF_Phase_Correction_Vec_1_draw(U,disPosu,tau,thre)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if nargin < 4, thre = 1e10; end
nDis = numel(disPosu);
cntd = 1;
U2 = U;
while (cntd<=nDis)
    if cntd<nDis
        edPoint = disPosu(cntd+1)-1;
    else
        edPoint = numel(U);
    end
    num = 0;
    if cntd > 1
        aug = round(((U(disPosu(cntd)-1)-U(disPosu(cntd)))+tau/2)/tau)*tau;
        vec = aug+(-2:2)*tau+U(disPosu(cntd));
        [~,pos] = min(abs(vec-U(disPosu(cntd)-1)));
        U(disPosu(cntd)) = vec(pos);
    end
    der1st = 0;
    mk2nd = 0;
    cnt = disPosu(cntd)+1;
    ct = 1;
    
    
    while (cnt<=edPoint & ct)
        pic = figure;plot(U2,'-sb','LineWidth',2);hold on;
        plot(U(1:cnt-1,:),'-sr','LineWidth',2);
        axis square;axis([0,32,-5,40]);
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/pic/p%d.eps',cnt);
        saveas(pic,tit,'epsc');
   
        
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
            if abs(der2nd)>thre
                if cntd<nDis
                    disPosu = [disPosu(1:cntd),cnt,disPosu(cntd+1:end)];
                else
                    disPosu = [disPosu(1:cntd),cnt];
                end
                edPoint = cnt - 1;
                ct = 0;
                nDis = nDis + 1;
            end
        else
            [~,pos] = min(abs(vec-U(cnt-1)-der1st-mk2nd));
            U(cnt) = vec(pos);
        end
        der1st = U(cnt)-U(cnt-1);
        cnt = cnt + 1;
    end
    cntd = cntd + 1;
end
end







