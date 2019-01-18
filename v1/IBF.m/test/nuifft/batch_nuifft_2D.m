log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

NGlist = [4 6 8];
tollist = [1e-4 1e-6 1e-8];

for N = 2.^(5:1:8)
    fid = fopen([log_path 'Factor_nuifft_2D_' num2str(N) '.log'],'a+');
    for iter = 1:length(NGlist)
        NG = NGlist(iter);
        tol = tollist(iter);
        run_nuifft_2D(N, NG, tol, fid);
        fprintf('Nuifft, N %4d, NG %2d finished.\n',N,NG);
    end
    fclose(fid);
end
