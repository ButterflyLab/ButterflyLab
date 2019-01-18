function M = mkmommat2(edges, verts, fintg,rd,cl)
% M = mkmommat2(edges, verts, fintg)
%
% Computes the moment matrix. The matrix size in NxN, where N is the number
% of the boundary elements. The entry M(i,j) corresponds to the integral of
% the Green's function evaluated over boundary element j tested against
% basis function i.
%  Params:
%    edges  - num_of_edges-by-2 matrix of the indices of the edge endpoint
%             vertices
%    verts  - num_of_verts-by-2 matrix of the coordinates of the vertices.
%    fintg  - handle of the function which evaluates inetgral of the
%             Green's function. The function should accept two
%             parameters: source and observation segments.
%
if nargin < 4
    Ne = size(edges,1);
    Nv = size(verts,1);
    rd = 1:Ne;
    cl = 1:Nv;
else
    Ne = numel(rd);
    Nv = numel(cl);
end
[m,n] = ndgrid(rd,cl);

m = m(:);
n = n(:);
rsrc = cat(3, verts(edges(n,1),:), verts(edges(n,2),:));
robs = cat(3, verts(edges(m,1),:), verts(edges(m,2),:));

M = zeros(Ne,Nv);
M(:) = fintg(rsrc,robs);
