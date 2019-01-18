function v = intg_lapdl2d(rsrc,robs)
% v = intg_lapdl2d(rsrc,robs)
%
% Evaluates integral of the double layer potential over the boundary
% segments, where the Green's function is the one associated with
% the 2D Laplace equation which is is given by
%   -1/(2*pi)*ln(r/c)
% The double layer potential is the normal derivative of the Green's
% function given by the following formula
%   g*(R) = -1/(2*pi*r)*dr/dR*n
% where r is distance, R is position vector, n is normal, dr/dR = R/r.
% Notice that this function identifies singular source-observation pairs
% (ones where source coincides with observation) by comparison of the
% segment endpoints which only allows to identify the cases where the
% source segment exactly matches the observation one, any other singular
% configurations are ignored and may cause undefined behavior.
%
% Params:
%  rsrc   - source segments, N-by-2-by-2 array, idx-XYZ-begin/end.
%  robs   - observation segments, N-by-2-by-2 array, idx-XYZ-begin/end.
%           

% Number of edges
N = size(rsrc,1);

% Edge vectors. N-by-2 array. First index is the edge index, second
% one is X-Y.
edges = rsrc(:,:,2)-rsrc(:,:,1);

% Edge lengths. Column vector of length N.
l = sqrt(sum(edges.^2,2));

% Edge tangentials - normalized edge vectors.
t = edges ./ l(:,ones(1,2));

% Outer normals. The outer polygon is CCW, the holes are CW, thus
% to get the outer normals we need to rotate the edge tangentials
% 90 degrees clockwise.
n = [ t(:,2) -t(:,1) ];

% Find singular pairs
sidx = find(all(all(abs(rsrc-robs)<1e-12, 3), 2));

% Edge centers, N-by-2
obsc = sum(robs, 3)*0.5;

r = rsrc(:,:,1) - obsc;
l0 = 1; % An arbitrary constant; makes the log argument dimensionless
l0_2 = l0*l0;

nr = dot(n,r,2);
absnr = abs(nr);

tr = dot(t,r,2);
signnr = sign(nr);
t1 = atan((l+tr)./absnr);
t2 = atan(tr./absnr);
v = signnr.*(t1-t2)*-1/(pi*2);

v(sidx) = 0;
