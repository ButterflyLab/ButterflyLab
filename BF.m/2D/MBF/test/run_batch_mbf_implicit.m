addpath('../src/');
log_path = './log/';
data_path = './data/';

if(~exist(data_path, 'dir'))
    mkdir(data_path);
end
if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

func_list = {'fun0'};
for func_i = 1:length(func_list)
    func_name = func_list{func_i};
    for N = 2.^(6:7)
        fid = fopen([log_path 'Factor_' func_name '_' func_name '_' num2str(N) '.log'],'a+');
        Factor1 = load([data_path 'Factor_' func_name '_' num2str(N) '_' num2str(28) '.mat'],'Factor');
        Factor2 = Factor1;
        fun = @(x)apply_mbf(Factor2.Factor,...
            reshape(fft2(reshape(apply_mbf(Factor1.Factor,x),N,N,[]))/N,N^2,[]));
        fun_adj = @(y)apply_mbf_adj(Factor1.Factor,...
            reshape(ifft2(reshape(apply_mbf_adj(Factor2.Factor,y),N,N,[]))*N,N^2,[]));
        for mR = [16 24]
            tol = 10^(-mR/8);
            run_mbf_implicit(N, fun, fun_adj, mR, tol, fid);
            fprintf('Func %s, N %4d, mR %2d finished.\n',func_list{func_i},N,mR);
        end
        fclose(fid);
    end
end
