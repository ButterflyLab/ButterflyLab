function opc = op_count(Factor)
% This code estimate the operation count according to the number of nonzero
% entries in the CUR butterfly factorization.
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

opc = nnz(Factor.V);

for i=length(Factor.HTol):-1:1
    opc = opc + nnz(Factor.HTol{i});
end

opc = opc + nnz(Factor.M);

for i=1:length(Factor.GTol)
    opc = opc + nnz(Factor.GTol{i});
end

opc = opc + nnz(Factor.U);

end
