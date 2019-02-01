function  FacOut = HSSBF_fwd_assemble(Factor)
% HSSBF_fwd_assemble assebles a HSSBF matrix stored in a data-sparse HSSBF
% format in Factor. Returns a structure storing matrix factors of the HSSBF
% matrix. Complexity: O(N\log(N)) in both operation and memory.
%
% Input:
% Factor: a structure storing the data-sparse representation of the matrix
%        A(x,k), and some auxiliary variables.
% Factor.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% Factor.sz: the size of A
% Factor.szsub: the size of A11
% Factor.A11
% Factor.A12
% Factor.A21
% Factor.A22
%
% Output:
% FacOut, a cell structure storing matrix factors of the HSSBF matrix, i.e.
%      A \aprox FacOut{1}*FacOut{2}*...*FacOut{2n+1}
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

if Factor.isLeaf
    A = [Factor.A11,Factor.A12;Factor.A21,Factor.A22];
    FacOut = cell(1,1);
    FacOut{1} = sparse(A);
else
    fac1 = HSSBF_fwd_assemble(Factor.A11);
    fac2 = HSSBF_fwd_assemble(Factor.A22);
    f1 = BF_rename(Factor.A21);
    f2 = BF_rename(Factor.A12);
    n1 = length(fac1);
    n2 = length(f1);
    if n1 == n2 % case 1
        n = n1;
        FacOut = cell(1,n);
        L1 = blkdiag(fac1{1},f1{1});
        R1 = blkdiag(f2{1},fac2{1});
        FacOut{1} = [L1,R1]; % permute first factor
        L2 = [fac1{n1};f1{length(f1)}];
        R2 = [f2{length(f2)};fac2{n1}];
        FacOut{n} = blkdiag(L2,R2);
        for cnt = 1:n-2
            FacOut{cnt+1} = blkdiag(fac1{cnt+1},f1{cnt+1},f2{cnt+1},fac2{cnt+1});
        end
    else % case 2
        n = n1+2;
        FacOut = cell(1,n);
        L1 = blkdiag(speye(Factor.szsub),f1{1});
        R1 = blkdiag(f2{1},speye(Factor.sz-Factor.szsub));
        FacOut{1} = [L1,R1]; % permute first factor
        L2 = [speye(Factor.szsub);f1{length(f1)}];
        R2 = [f2{length(f2)};speye(Factor.sz-Factor.szsub)];
        FacOut{n} = blkdiag(L2,R2);
        for cnt = 1:n-2
            FacOut{cnt+1} = blkdiag(fac1{cnt},f1{cnt+1},f2{cnt+1},fac2{cnt});
        end
    end
end