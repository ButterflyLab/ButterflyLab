function [grid,Lgrid] = BF_grid(x,npx,gridIn,box)
% This code arranges the grid points for the butterfly factorization.
%
% x - the x-th box
% npx - the total number of boxes
% gridIn - the unscaled grid on [0,1]
% box - [box(1),box(2)] is the new domain of the scaled grid points
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

Dim = size(x,2);
len = (box(2,:)-box(1,:))./npx;
offset = box(1,:)+(x-1).*len;
Lgrid = gridIn*len+ones(size(gridIn))*offset;

I=cell(Dim,1);
for di=1:Dim
    I{di}=Lgrid(:,di);
end
[I{1:Dim}]=ndgrid(I{:});

grid = zeros(length(gridIn)^Dim,Dim);
for di=1:Dim
    grid(:,di)=I{di}(:);
end

end