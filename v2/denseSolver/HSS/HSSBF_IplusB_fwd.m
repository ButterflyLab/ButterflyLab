function  Factor = HSSBF_IplusB_fwd(Bfactor,lsz)
% HSSBF_IplusB_fwd compresses a HSSBF matrix A=(I+B), where B is stored in
% a data-sparse BF format in Bfactor.
%
% Input:
% Bfactor: a data-sparse BF format of the complementary low-rank matrix B.
% lsz:  a size parameter; when the matrix A has a size <= lsz, we
%       don't compress
%
% Output:
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
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

M = size(Bfactor.U,1); N = size(Bfactor.V,2);
Factor = struct('isLeaf',[],'sz',[],'A11',[],'A12',[],'A21',[],'A22',[]);
if M~=N
    error('matrix is not square');
end
if M>lsz
    Factor.sz = M;
    [B11,B12,B21,B22,flag] = BF_divide1D(Bfactor);
    if flag == 0
        Factor.szsub = size(B11.U,1);
        Factor.isLeaf = 0;
        Factor.A11 = HSSBF_IplusB_fwd(B11,lsz);
        Factor.A22 = HSSBF_IplusB_fwd(B22,lsz);
        Factor.A12 = B12;
        Factor.A21 = B21;
    else
        Factor.szsub = size(B11,1);
        Factor.isLeaf = 1;
        Factor.A11 = B11+eye(Factor.szsub);
        Factor.A22 = B22+eye(M-Factor.szsub);
        Factor.A12 = B12;
        Factor.A21 = B21;
    end
else
    Factor.sz = M;
    [B11,B12,B21,B22,flag] = BF_divide1D(Bfactor);
    if flag == 0
        Factor.szsub = size(B11.U,1);
        Factor.isLeaf = 1;
        Factor.A11 = BF_sp2den(B11)+eye(Factor.szsub);
        Factor.A22 = BF_sp2den(B22)+eye(M-Factor.szsub);
        Factor.A12 = BF_sp2den(B12);
        Factor.A21 = BF_sp2den(B21);
    else
        Factor.szsub = size(B11,1);
        Factor.isLeaf = 1;
        Factor.A11 = B11+eye(Factor.szsub);
        Factor.A22 = B22+eye(M-Factor.szsub);
        Factor.A12 = B12;
        Factor.A21 = B21;
    end
end