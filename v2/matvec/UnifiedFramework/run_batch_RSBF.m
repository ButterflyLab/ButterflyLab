
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

func_list = {'fun1','fun2','funH'};
method = 1;
mR = 20;
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(10:2:18)
        fid = fopen(['./results/RSBF_' func_name '_' num2str(N) '.log'],'a+');
        for NG = 6:2:12
            tol = 1e-12;
            run_RSBF(N, func_name, NG, mR, tol, fid, method);
            fprintf('Func %s, N %4d, NG %2d, mR %2d finished.\n',func_list{func_i},N,NG,mR);
        end
        fclose(fid);
    end
end
