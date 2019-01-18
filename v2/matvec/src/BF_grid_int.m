function grid = BF_grid_int(gridIn,x,NG,npx,box,disSet)
% This code arranges the grid points for the butterfly factorization.
%
% gridIn - the unscaled Chebyshev grid on [0,1]
% x - the x-th box
% NG - the number of points in a grid
% npx - the total number of boxes
% box - [box(1),box(2)] is the new domain of the scaled grid points
% disSet - a set of indices that need to exclude due to the fact that the
% phase function is not smooth at these locations.
%
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if nargin < 6
    numDis = 0;
else
    numDis = numel(disSet);
end

len = (box(2)-box(1)+1)./npx;
st = ceil(box(1)+(x-1).*len);
ed = floor(box(1)+x.*len-1);
sz = numel(gridIn);
grid = round(gridIn*((ed-st+1)-sz) + (0:sz-1)') + st;
grid = unique(grid);
if numDis > 0, setdiff(grid,disSet); end
end