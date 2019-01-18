log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

%WarmUP
ifft(rand(128,1));

NGlist = [4 6 8];
tollist = [1e-4 1e-7 1e-10];

for N = 2.^(8:2:18)
    x = rand(N,1);
    tic;
    ifx = ifft(x);
    Tifft = toc;
    fprintf('IFFT, N %4d cost %.2e secs.\n',N,Tifft);
    fid = fopen([log_path 'Factor_nuifft_1D_' num2str(N) '.log'],'a+');
    for iter = 1:length(NGlist)
        NG = NGlist(iter);
        tol = tollist(iter);
        run_nuifft_1D(N, NG, tol, fid);
        fprintf('Nuifft, N %4d, NG %2d finished.\n',N,NG);
    end
    fclose(fid);
end
