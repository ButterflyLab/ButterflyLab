log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'funF'};
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(6:9)
        fid = fopen([log_path 'Factor_' func_name '_' num2str(N) '.log'],'a+');
        for mR = 6:2:10
            tol = 1e-12;
            run_bf_explicit(N, func_name, mR, tol, fid, mR==6^2);
            fprintf('Func %s, N %4d, mR %2d finished.\n',func_list{func_i},N,mR);
        end
        fclose(fid);
    end
end
