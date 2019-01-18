function [fun,xx,kk,xbox,kbox,npx,npk] = BF_prepbox(fun_org,xx,kk,tag)
% This code arrange the domain of the kernel function or matrix for the
% butterfly factorization.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

Dim = size(xx,2);
xbox = zeros(2,Dim);
kbox = zeros(2,Dim);

if isscalar(tag)
    N = tag^(1/Dim)*ones(1,Dim);
    [Nx,~] = size(xx);
    [Nk,~] = size(kk);
    Nx = Nx^(1/Dim)*ones(1,Dim);
    Nk = Nk^(1/Dim)*ones(1,Dim);
    npx = 2.^ceil(log2(sqrt(N)));
    npk = 2.^floor(log2(sqrt(N)));
    xbox(1,:) = min(xx);
    xbox(2,:) = min(xx)+(max(xx)-min(xx))./Nx.*N;
    kbox(1,:) = min(kk);
    kbox(2,:) = min(kk)+(max(kk)-min(kk))./Nk.*N;
    fun = @(x,k)fun_org(x,k);
    return;
end

if strcmpi(tag,'regular')
    [Nx,~] = size(xx);
    [Nk,~] = size(kk);
    Nx = Nx^(1/Dim)*ones(1,Dim);
    Nk = Nk^(1/Dim)*ones(1,Dim);
    npx = 2.^ceil(log2(sqrt(Nx)));
    npk = 2.^floor(log2(sqrt(Nk)));
    xbox(1,:) = min(xx);
    xbox(2,:) = max(xx)+(max(xx)-min(xx))./Nx;
    kbox(1,:) = min(kk);
    kbox(2,:) = max(kk)+(max(kk)-min(kk))./Nk;
    fun = @(x,k)fun_org(x,k);
    return;
end

if strcmpi(tag,'polar')
    [Nx,~  ] = size(xx);
    kk = BF_k2p(kk);
    [Nk,~  ] = size(kk);
    Nx = Nx^(1/Dim)*ones(1,Dim);
    Nk = Nk^(1/Dim)*ones(1,Dim);
    npx = 2.^ceil(log2(sqrt(Nx)));
    npk = 2.^floor(log2(sqrt(Nk)));
    npk(2) = 4*npk(2);
    if Dim == 3
        npk(3) = 2*npk(3);
    end
    xbox(1,:) = min(xx);
    xbox(2,:) = max(xx)+(max(xx)-min(xx))./Nx;
    kbox(1,:) = min(kk);
    kbox(2,:) = max(kk)+(max(kk)-min(kk))./Nk;
    fun = @(x,p)fun_org(x,BF_p2k(p));
    return;
end

if strcmpi(tag,'multiscale')
    [Nx,~  ] = size(xx);
    [Nk,~  ] = size(kk);
    Nk = Nk/3*4;
    Nx = Nx^(1/Dim)*ones(1,Dim);
    Nk = Nk^(1/Dim)*ones(1,Dim);
    npx = 2.^ceil(log2(sqrt(Nk)));
    npk = 2.^floor(log2(sqrt(Nk)));
    xbox(1,:) = min(xx);
    xbox(2,:) = max(xx)+(max(xx)-min(xx))./Nx;
    kbox(1,:) = min(kk);
    kbox(2,:) = max(kk)+(max(kk)-min(kk))./Nk;
    fun = @(x,k)fun_org(x,k);
    return;
end

end