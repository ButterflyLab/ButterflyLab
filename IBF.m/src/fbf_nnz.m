function nnz = fbf_nnz(levels,Mcell,GTolcell,Ucell,HTolcell,Vcell)

nnz = 0.0;

%---------------------------------------------------------------
%   M nnz
for it = 1:length(Mcell)
    nnz = nnz + numel(Mcell{it});
end


%---------------------------------------------------------------
%   G nnz
for ell = 1:levels
    for it = 1:length(GTolcell{ell})
        nnz = nnz + numel(GTolcell{ell}{it});
    end
end


%---------------------------------------------------------------
%   U nnz
for it = 1:length(Ucell)
    nnz = nnz + numel(Ucell{it});
end


%---------------------------------------------------------------
%   H nnz
for ell = 1:levels
    for it = 1:length(HTolcell{ell})
        nnz = nnz + numel(HTolcell{ell}{it});
    end
end


%---------------------------------------------------------------
%   V nnz
for it = 1:length(Vcell)
    nnz = nnz + numel(Vcell{it});
end

end