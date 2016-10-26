log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

NGlist = [2 8];
tollist = [1e-1 1e-9];

for N = 2.^(9)
    fid = fopen([log_path 'Factor_nufft_2D_' num2str(N) '.log'],'a+');
    for iter = 1:length(NGlist)
        NG = NGlist(iter);
        tol = tollist(iter);
        run_nufft_2D(N, NG, tol, fid);
        fprintf('Nufft, N %4d, NG %2d finished.\n',N,NG);
    end
    fclose(fid);
end
