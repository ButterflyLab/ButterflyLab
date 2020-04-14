% Multidimensional Phase Recovery and Interpolative Decomposition 
% Butterfly Factorization, preprint, 2019.
% Copyright 2018 by Ze Chen, Juan Zhang, Ken Ho, and Haizhao Yang

% Z-order curve for the partition of the problem domain.

function order = BF_Zorder(N, dim)

switch dim
    case 1
        order = (1:N)';
    case 2
        bits = log2(N);
        bin = dec2bin(0:N-1,bits);
        order = sub2ind(2.^([bits, bits]/2), ...
            bin2dec(bin(:,2:2:bits))+1, ...
            bin2dec(bin(:,1:2:bits-1))+1);
    case 3
        bits = log2(N);
        bin = dec2bin(0:N-1,bits);
        order = sub2ind(2.^([bits, bits, bits]/3), ...
            bin2dec(bin(:,3:3:bits))+1, ...
            bin2dec(bin(:,2:3:bits-1))+1, ...
            bin2dec(bin(:,1:3:bits-2))+1);
end