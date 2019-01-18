function opc = op_count(Factor)

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
