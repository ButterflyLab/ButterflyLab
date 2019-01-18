function [xbox,kbox,npx,npk] = BF_prepbox_int(Nx,Nk,tag)
% This code arrange the domain of the kernel function or matrix for the
% butterfly factorization.
%
% Copyright 2017 by Haizhao Yang

if strcmpi(tag,'regular')
    xbox = zeros(2,1);
    kbox = zeros(2,1);
    npx = 2.^ceil(log2(sqrt(Nx)));
    npk = 2.^floor(log2(sqrt(Nk)));
    xbox(1,:) = 1;
    xbox(2,:) = Nx;
    kbox(1,:) = 1;
    kbox(2,:) = Nk;
    return;
end

if strcmpi(tag,'multiscale')
    xbox = zeros(2,1);
    kbox = zeros(2,1);
    %Nk = Nk/3*4;
    npx = 2.^ceil(log2(sqrt(Nk)));
    npk = 2.^floor(log2(sqrt(Nk)));
    xbox(1,:) = 1;
    xbox(2,:) = Nx;
    kbox(1,:) = 1;
    kbox(2,:) = Nk;
    return;
end
end