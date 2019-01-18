function mask = BF_nnz_position(BF)
% This code identify the nonzero mask of the butterfly factorization.
%
% Input:
% BF - the butterfly factorization storing sparse matrices U, GTol{},M,
%      HTol{}, and V such that
%     BF.U*BF.HTol{Gl}*...*BF.HTol{1}*BF.M*BF.HTol{1}*...*BF.HTol{Hl}*BF.V
%     is the butterfly factorization of a complementary matrix.
%
% Output:
% mask - a cell of size 4 by 3+Hl+Gl, where the (1,j)-th cell stores the
%        row indices of nonzero entries in the j-th factor of BF, and the (2,j)-th
%        cell stores the column indices of nonzero entries in the j-th factor
%         of BF, the j-th factor is of size mask{3,j} by mask{4,j}.
%
% Copyright 2018 by Haizhao Yang

L = 3+length(BF.HTol)+length(BF.GTol);
mask = cell(2,L);
[mask{1,L},mask{2,L},~] = find(BF.V);
[mask{3,L},mask{4,L}] = size(BF.V);

Hl = length(BF.HTol);
for i=Hl:-1:1
    [mask{1,L-(Hl-i)-1},mask{2,L-(Hl-i)-1},~] = find(BF.HTol{i});
    [mask{3,L-(Hl-i)-1},mask{4,L-(Hl-i)-1}] = size(BF.HTol{i});
end

[mask{1,L-1-Hl},mask{2,L-1-Hl},~] = find(BF.M);
[mask{3,L-1-Hl},mask{4,L-1-Hl}] = size(BF.M);

Gl = length(BF.GTol);
for i=1:Gl
    [mask{1,L-1-Hl-i},mask{2,L-1-Hl-i},~] = find(BF.GTol{i});
    [mask{3,L-1-Hl-i},mask{4,L-1-Hl-i}] = size(BF.GTol{i});
end

[mask{1,1},mask{2,1},~] = find(BF.U);
[mask{3,1},mask{4,1}] = size(BF.U);

end
