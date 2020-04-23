clear; clc;
close all hidden;
rng('shuffle');

N = 2 .^ (7 : 17);
fun_m = @(n) round(1.0*n);

% Set up parameters
tol = 1e-10;
r = 100;
mR = 30;
nn = 256;
datatry = 10;

tstart = tic;
startT = char(datetime('now'));
fid = fopen(['./results/' startT '.log'], 'a+');
fprintf(fid, ['\n\n Starting from ', startT, '\n']);
fprintf(fid, [' epsilon: ', num2str(tol), '\tk_epsilon: ', num2str(r), ...
    '\tr_epsilon: ', num2str(mR), '\tmin block: ', num2str(nn), '\n']);
fprintf(fid, ['N & m & Fwd Err & Inv Err & Total Err & ' ...
    'Fwd Mat/Fac & Fwd Dir/App & Inv Mat/Fac & Inv Dir/App\n']);

cm = length(N);
Err = zeros(cm, datatry); Err_f = Err; Err_i = Err;
AppT_f = Err; AppT_i = Err;
FacT_f = zeros(cm, 1); FacT_i = FacT_f;
MatT_f = FacT_f; MatT_i = FacT_f;
DirT_f = FacT_f; DirT_i = FacT_f;
pf = cell(cm, 1);

pid = fopen(['./results/processing_' startT '.log'], 'a+');
for c = 1 : cm
    
    profile on;
    
    n = N(c);
    m = fun_m(n);
    
    fprintf(pid, '\nN: %i, test m: %i, matrix size: %i*%i\n', n, m, 2*n, 2*n-m);
    [theta, weight] = legenQuad(2*n);
    theta = acos(real(theta(:))); weight = real(weight(:));
    fun = @(x,k) ALegendrefun1(x,k,m,theta,(m:2*n-1)',weight);
    
    proc = (c-1)/cm;
    fprintf(pid, ['%2.2f%%\t%2.0f/%2.0f\t', char(datetime('now')),'\n'], ...
        proc*100, c, cm);
    
    a = randn(2*n-abs(m), datatry);
    
    % a2g
    [g, FacT_f(c), AppT_f(c,:)] = ...
        BF_ALegendre_fwd(abs(m),2*n-1,a,theta,weight,r,mR,tol,nn,datatry);
    AppT_f(c,:) = AppT_f(c,:) / datatry;
    
    cut = find(theta(n+1:n*2) > real(asin(sqrt(m^2-1/4)/...
        ((2*n-1)+1/2))),1,'last') + n;
    NC = min(128, cut - n);
    ind = sort(randperm(cut - n, NC) + n);
    t1 = tic;
    KF = zeros(NC*2, 2*n-m);
    KF(end/2+1:end, :) = fun(ind, 1:(2*n-m));
    KF(1:end/2, 1:2:end) = KF(end:-1:end/2+1, 1:2:end);
    KF(1:end/2, 2:2:end) = -KF(end:-1:end/2+1, 2:2:end);
    MatT_f(c) = toc(t1);
    MatT_f(c) = MatT_f(c) * n / NC;
    t2 = tic;
    gext = KF*a;
    DirT_f(c) = toc(t2);
    DirT_f(c) = DirT_f(c) * n / NC / datatry;
    Err_f(c,:) = vecnorm(gext - g([flip(2*n+1-ind), ind],:)) ./ vecnorm(gext);
    
    % g2a
    [aa, FacT_i(c), AppT_i(c,:)] = ...
        BF_ALegendre_inv(abs(m),2*n-1,g,theta,weight,r,mR,tol,nn,datatry);
    AppT_i(c,:) = AppT_i(c,:) / datatry;
    
    NC = min(256, 2*n-m);
    ind = sort(randperm(2*n-m, NC));
    iodd = find(mod(ind,2) == 1); ieven = find(mod(ind,2) == 0);
    t1 = tic;
    KI = zeros(NC, 2*n);
    KI(:, end/2+1:end) = fun(n+1:2*n, ind)';
    KI(iodd, 1:end/2) = KI(iodd, end:-1:end/2+1);
    KI(ieven, 1:end/2) = -KI(ieven, end:-1:end/2+1);
    MatT_i(c) = toc(t1);
    MatT_i(c) = MatT_i(c) * (2*n-m) / NC;
    t2 = tic;
    aaext = KI*g;
    DirT_i(c) = toc(t2);
    DirT_i(c) = DirT_i(c) * (2*n-m) / NC / datatry;
    Err_i(c,:) = vecnorm(aaext - aa(ind,:)) ./ vecnorm(aaext);
    
    Err(c,:) = vecnorm(a - aa) ./ vecnorm(a);
    
    profile off;
    pf{c} = profile('info');
    
    save(['./results/' startT '.mat'],'N','Err','Err_f','Err_i',...
        'FacT_f','FacT_i','AppT_f','AppT_i','MatT_f','MatT_i',...
        'DirT_f','DirT_i','pf');
    
    LT_scaling(startT, 3, 5, 0)
    
    fprintf(fid, ['%i & %i ' repmat('& %.2e ',1,7) '\\\\\n'], ...
        n, m, BF_trimdata(Err_f(c,:)), BF_trimdata(Err_i(c,:)), ...
        BF_trimdata(Err(c,:)), MatT_f(c)/FacT_f(c), ...
        DirT_f(c)/BF_trimdata(AppT_f(c,:)), ...
        MatT_i(c)/FacT_i(c), ...
        DirT_i(c)/BF_trimdata(AppT_i(c,:)));
    
end

T = toc(tstart) / 60;
fprintf(fid, ' Running time = %.2f min\n', T);
fprintf(pid, '\nFinished. Running time = %.2f min\n', T);
