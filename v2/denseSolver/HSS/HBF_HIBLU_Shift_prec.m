function x = HBF_HIBLU_Shift_prec(Zfun,b,restart,tolsol,maxit)
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

[x,~,err,iter,~] = gmres(Zfun,b,restart,tolsol,maxit);
(iter(1)-1)*restart+iter(2)
err
end