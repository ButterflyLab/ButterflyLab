function f = BF_interpolation2D(V,xgd,ygd,x,y)
% 2D Lagrange interpolation
%
% Input:
% xgd - grid points in x
% ygd - grid points in y
% V - function values on the grid
% x - new points in x of interest
% y - new points in y of interest
%
% Output
% function values f(x,y)

nx = numel(x); ny = numel(y);
f = zeros(nx,ny);
[Nx,Ny] = size(V);
ellx = zeros(Nx,1);
elly = zeros(Ny,1);
for cntx = 1:nx
    xcur = x(cntx);
    % interpolation in x
    for ith = 1:Nx
        ellx(ith) = lagrange_basis_function_1d (Nx-1,xgd,ith,xcur);
    end
    Vx = ellx.'*V;
    for cnty = 1:ny
        ycur = y(cnty);
        % interpolation in y
        for ith = 1:Ny
            elly(ith) = lagrange_basis_function_1d (Ny-1,ygd,ith,ycur);
        end
        f(cntx,cnty) = Vx*elly;
    end
end
