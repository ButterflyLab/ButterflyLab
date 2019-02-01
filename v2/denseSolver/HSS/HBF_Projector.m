function P0 = HBF_Projector(Z,N,numPoles,R,numEigVal,thre)
% Use contour integral to generate a projection P corresponding to the
% eigensubspace of Z with eigenvalues inside a circle centered at the
% origin with a radius R.
%
% Z, a function handle such that Z(z,b) solves (Z-z*I)x=b with the HIBLU
% preconditioner
% N? size of Z
% numPoles, the number of poles to discretize the contour
% R, the radius of the contour
% P0, the spectral projector corresponding to the eigenvalues inside the
% given circle
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

%% use the trapizoidal rule to discretize the contour
Emin = -R; Emax = R;
[lgx,lgw] = HBF_lgwt(numPoles,-1,1);
rad = (Emax-Emin)/2; theta = -(pi/2)*(lgx-1); %theta in (0,pi)
ex = exp(i*theta);
z = (Emax+Emin)/2+rad*ex; % positions of the poles
w = rad*lgw.*ex/2;  % weights of the poles


%% solve h_j=(H-z_j)^{-1}b for the quadrature points z_j using HIBLU
hh = cell(1,numPoles);
% use exact method to solve linear systems
for cnt = 1:numPoles
    hh{cnt} = @(b) Z(z(cnt),b);
end

%----------------------------------------------------------------------
% compute an approximated eigenspace corresponding to numEigVal
% smallest eigenvalues of H using the contour integral

b = rand(N,numEigVal);
Pb = zeros(size(b));
for cnt=1:numel(w)
    Pb = Pb+real(w(cnt)*hh{1,cnt}(b));
end

[Pb,R] = qr(Pb,0);
pos = find(diag(R)>thre);
Pb = Pb(:,pos);
P0 = @(x) Pb*(Pb'*x);
end

