log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'fun0'};
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(6:8)
        fid = fopen([log_path 'Factor_' func_name '_' num2str(N) '.log'],'a+');
        for mR = [6 14 22]
            tol = 10^(-(mR+2)/8-1);
            run_pbf_explicit(N, func_name, mR, tol, fid, mR==22);
            fprintf('Func %s, N %4d, mR %2d finished.\n',func_list{func_i},N,mR);
        end
        fclose(fid);
    end
end
