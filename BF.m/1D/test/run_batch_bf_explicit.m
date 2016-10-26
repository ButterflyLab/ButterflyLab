log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'fun0'};
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(8:2:18)
        fid = fopen([log_path 'Factor_' func_name '_' num2str(N) '.log'],'a+');
        for mR = 4:2:8
            tol = 1e-12;
            run_bf_explicit(N, func_name, mR, tol, fid);
            fprintf('Func %s, N %4d, mR %2d finished.\n',func_list{func_i},N,mR);
        end
        fclose(fid);
    end
end
