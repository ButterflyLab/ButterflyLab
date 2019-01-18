function tol = BF_nnz_show(factor)
% This code estimate the number of nonzero entries in the CUR butterfly
% factorization for pre-allocation.
%
% Copyright 2018 Haizhao Yang


switch numel(fieldnames(factor))
    case 5
        tol = 0;
        t = nnz(factor.M)
        tol = tol + t;
        for ell = 1:length(factor.GTol)
            t = nnz(factor.GTol{ell})
            tol = tol + t;
        end
        t = nnz(factor.U)
        tol = tol + t;
        for ell = 1:length(factor.HTol)
            t = nnz(factor.HTol{ell})
            tol = tol + t;
        end
        t = nnz(factor.V)
        tol = tol + t;
    case 4
        tol = 0;
        r = length(factor.U);
        for ii = 1:r
            tol = tol + nnz(factor.V{ii});
        end
        tol = tol + nnz(factor.S);
        for ii = r:-1:1
            tol = tol + nnz(factor.U{ii});
        end
    case 3
        tol = 0;
        r = length(factor.U);
        for ii = 1:r
            tol = tol + nnz(factor.V{ii});
        end
        tol = tol + nnz(factor.S);
        for ii = r:-1:1
            tol = tol + nnz(factor.U{ii});
        end
end
end