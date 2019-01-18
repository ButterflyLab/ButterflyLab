
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

func_list = {'fun2'};
mR = 20;
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(10:2:18)
        fid = fopen(['./results/MVBFCUR_' func_name '_' num2str(N) '.log'],'a+');
        for NG = 6:2:12
            load(['Factor_' func_name '_' num2str(N) '_' num2str(12) '_' num2str(mR) '.mat']);
            fun = @(x) BF_apply(Factor,fft(BF_apply(Factor,x))/sqrt(N));
            fun_adj = @(y) BF_adj_apply(Factor,ifft(BF_adj_apply(Factor,y))*sqrt(N));
            tol = 1e-12;
            run_MVBF(N, fun, fun_adj, NG, mR, tol, fid,2);
            fprintf('Func %s, N %4d, NG %2d, mR %2d finished.\n',func_list{func_i},N,NG,mR);
        end
        fclose(fid);
    end
end
