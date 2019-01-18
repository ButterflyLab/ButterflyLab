log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

NGlist = [2 6 9];
tollist = [5e-1 1e-4 1e-7];

func_list = {'fun0'};
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(5:1:8)
        fid = fopen([log_path 'Factor_' func_name '_2D_PBF_' num2str(N) '.log'],'a+');
        for iter = 1:length(NGlist)
            NG = NGlist(iter);
            tol = tollist(iter);
            run_fastpbf_2D(N, func_name, NG, tol, fid);
            fprintf('Func %s, N %4d, NG %2d finished.\n',func_list{func_i},N,NG);
        end
        fclose(fid);
    end
end
