function c = HSSBF_LUpred(Zu)

%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

Zl = tril(Zu);
Zr = triu(Zu);
Zr = Zr - diag(diag(Zr)-ones(size(Zu,1),1));
Zu = Zl\(Zu/Zr);
c = cond(Zu);
end