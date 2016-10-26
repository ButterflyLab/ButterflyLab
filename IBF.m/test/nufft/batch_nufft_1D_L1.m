log_path = './log/';

if(~exist(log_path, 'dir'))
    mkdir(log_path);
end

%WarmUP
fft(rand(128,1));

NGlist = [2 6];
tollist = [1e-1 1e-4];

for N = 2.^(20)
    x = rand(N,1);
    tic;
    fx = fft(x);
    Tfft = toc;
    fprintf('FFT, N %4d cost %.2e secs.\n',N,Tfft);
    fid = fopen([log_path 'Factor_nufft_1D_' num2str(N) '.log'],'a+');
    for iter = 1:length(NGlist)
        NG = NGlist(iter);
        tol = tollist(iter);
        run_nufft_1D(N, NG, tol, fid);
        fprintf('Nufft, N %4d, NG %2d finished.\n',N,NG);
    end
    fclose(fid);
end
