function  [L,U] = LUBF(Afac,x,k,rk,tol,mr,lsz,method,pha0)
% LUBF computes the hierarchical LU factorization of a HSSBF markix stored
% in Afac as a data-sparse HSSBF format. Complexity: O(N\log(N)) in
% both operation and memory.
%
% Input:
% Afac: a markix stored in a HSSBF format.
% x and k: denote the row and column indices of the submarkix of A to be
%       factorizaed, i.e. LU factorize the submarkix A(x,k).
% rk:    the maximum rank of the low-rank approximation
% tol:  the accuracy parameter of the low-rank approximation. If we can
%       achieve a low-rank approximation with a rank smaller than rk and an
%       accuracy tol, we will use a smaller rank.
% mr:   the maximum rank of the amplitude and phase recovery
% lsz:  a size parameter; when the markix A has a size <= lsz, we
%       don't compress
% method: specify which BF to computer HSSBF
%         1: IDBF
%         2: CURBF
%         3: CURBF2
% pha0: phase function in the whole domain
%
% Output:
% L and U are in a HSSBF format
% L: a srkucture storing the data-sparse representation of the markix L
% in [L,U] = lu(A), and some auxiliary variables.
% L.isLeaf: 1, if all markices are dense
%           0, if all markices are data-sparse
% L.sz: the size of L
% L.szsub: the size of A11
% L.A11
% L.A12
% L.A21
% L.A22
%
% U: a srkucture storing the data-sparse representation of the markix U
% in [L,U] = lu(A), and some auxiliary variables.
% U.isLeaf: 1, if all markices are dense
%           0, if all markices are data-sparse
% U.sz: the size of U
% U.szsub: the size of A11
% U.A11
% U.A12
% U.A21
% U.A22
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

if nargin < 8, method = 3; end
if nargin < 9, pha0 = []; end

M = numel(x); N = numel(k);
L = struct('isLeaf',[],'sz',[],'szsub',[],'A11',[],'A12',[],'A21',[],'A22',[]);
U = struct('isLeaf',[],'sz',[],'szsub',[],'A11',[],'A12',[],'A21',[],'A22',[]);
if M~=N
    error('markix is not square');
end
if M>lsz
    L.sz = M; U.sz = N;
    L.szsub = numel(x(1:end/2)); U.szsub = numel(k(1:end/2));
    L.isLeaf = 0; U.isLeaf = 0;
    L.A12 = []; U.A21 = [];
    % step 1
    if isempty(pha0)
        [L.A11,U.A11] = LUBF(Afac.A11,x(1:end/2),k(1:end/2),rk,tol,mr,lsz,method);
    else
        pha = @(i,j) pha0(i,j);
        [L.A11,U.A11] = LUBF(Afac.A11,x(1:end/2),k(1:end/2),rk,tol,mr,lsz,method,pha);
    end
    % step 3
    % fast matvec for A21*inv(U11) and its conjugate rkanspose
    A21U11 = @(f) BF_apply(Afac.A21,LUBF_sol(U.A11, f,'U'));
    A21U11t = @(f) LUBF_adj_sol(U.A11, BF_adj_apply(Afac.A21,f),'U');
    N1 = L.sz-L.szsub; N2 = L.szsub;
    dim = 1;
    % low-rank approximation of the phase function
    if isempty(pha0)
        [Up,Sp,Vp] = MVBF_Lowrank_Phase(N1,N2,A21U11,A21U11t,tol,5*mr,mr,dim,pi,1);
    else
        pha = @(i,j) pha0(i+M/2,j);
        [Up,Sp,Vp] = MVBF_Lowrank_Phase_initial(N1,N2,A21U11,A21U11t,tol,5*mr,mr,dim,pi,pha,1);
    end
    % low-rank approximation of the amplitude function
    [Ua,Sa,Va] = MVBF_Lowrank_Amp(N1,N2,A21U11,A21U11t,tol,5*mr,mr,dim);
    % BF of L12
    Vp = Sp*Vp'; Va = Sa*Va';
    funk = @(t,s) (Ua(t,:)*Va(:,s)).*exp(1i*Up(t,:)*Vp(:,s));
    tt = 1:N1; tt = tt(:);
    ss = 1:N2; ss = ss(:);
    switch method
        case 1
            L.A21 = BF_IDBF(funk,tt,ss,max(rk,8),rk*5,tol);
        case 2
            L.A21 = CURBF(funk,tt,ss,rk,tol);
        case 3
            L.A21 = CURBF2(funk,tt,ss,max(rk,8),rk,tol);
    end
    
    
    % step 2
    % fast matvec for inv(L11)A12 and its conjugate rkanspose
    L11A12 = @(f) LUBF_sol(L.A11, BF_apply(Afac.A12,f),'L');
    L11A12t = @(f) BF_adj_apply(Afac.A12,LUBF_adj_sol(L.A11, f,'L'));
    N1 = L.szsub; N2 = L.sz-L.szsub;
    dim = 1;
    % low-rank approximation of the phase function
    if isempty(pha0)
        [Up,Sp,Vp] = MVBF_Lowrank_Phase(N1,N2,L11A12,L11A12t,tol,5*mr,mr,dim,pi,1);
    else
        pha = @(i,j) pha0(i,j+N/2);
        [Up,Sp,Vp] = MVBF_Lowrank_Phase_initial(N1,N2,L11A12,L11A12t,tol,5*mr,mr,dim,pi,pha,1);
    end
    % low-rank approximation of the amplitude function
    [Ua,Sa,Va] = MVBF_Lowrank_Amp(N1,N2,L11A12,L11A12t,tol,5*mr,mr,dim);
    % BF of U12
    Vp = Sp*Vp'; Va = Sa*Va';
    funk = @(t,s) (Ua(t,:)*Va(:,s)).*exp(1i*Up(t,:)*Vp(:,s));
    tt = 1:N1; tt = tt(:);
    ss = 1:N2; ss = ss(:);
    switch method
        case 1
            U.A12 = BF_IDBF(funk,tt,ss,max(rk,8),rk*5,tol);
        case 2
            U.A12 = CURBF(funk,tt,ss,rk,tol);
        case 3
            U.A12 = CURBF2(funk,tt,ss,max(rk,8),rk,tol);
    end
    
    % step 4
    % recover L21*U12
    
    L21U12 = @(f) BF_apply(L.A21, BF_apply(U.A12,f));
    L21U12t = @(f) BF_adj_apply(U.A12,BF_adj_apply(L.A21, f));
    N1 = L.sz-L.szsub; N2 = L.sz-L.szsub;
    dim = 1;
    % low-rank approximation of the phase function
    if isempty(pha0)
        [Up,Sp,Vp] = MVBF_Lowrank_Phase(N1,N2,L21U12,L21U12t,tol,5*mr,mr,dim,pi,1);
    else
        pha = @(i,j) pha0(i+M/2,j+N/2);
        [Up,Sp,Vp] = MVBF_Lowrank_Phase_initial(N1,N2,L21U12,L21U12t,tol,5*mr,mr,dim,pi,pha,1);
    end
    % low-rank approximation of the amplitude function
    [Ua,Sa,Va] = MVBF_Lowrank_Amp(N1,N2,L21U12,L21U12t,tol,5*mr,mr,dim);
    % BF of U12
    Vp = Sp*Vp'; Va = Sa*Va';
    % update A22
    if isempty(pha0)
        Afac.A22 = HSSBF_update(Afac.A22,Ua,Va,Up,Vp,x((end/2+1):end),k((end/2+1):end),rk,tol,mr,lsz,method);
        [L.A22,U.A22] = LUBF(Afac.A22,x((end/2+1):end),k((end/2+1):end),rk,tol,mr,lsz,method);
    else
        Afac.A22 = HSSBF_update(Afac.A22,Ua,Va,Up,Vp,x((end/2+1):end),k((end/2+1):end),rk,tol,mr,lsz,method,pha);
        [L.A22,U.A22] = LUBF(Afac.A22,x((end/2+1):end),k((end/2+1):end),rk,tol,mr,lsz,method,pha);
    end
else
    A = [Afac.A11,Afac.A12;Afac.A21,Afac.A22];
    [LL,UU] = lu(A);
    L.sz = M;
    L.szsub = numel(x(1:end/2));
    L.isLeaf = 1;
    L.A11 = LL(1:L.szsub,1:L.szsub);
    L.A22 = LL((L.szsub+1):end,(L.szsub+1):end);
    L.A12 = [];
    L.A21 = LL((L.szsub+1):end,1:L.szsub);
    U.sz = N;
    U.szsub = numel(k(1:end/2));
    U.isLeaf = 1;
    U.A11 = UU(1:L.szsub,1:L.szsub);
    U.A22 = UU((L.szsub+1):end,(L.szsub+1):end);
    U.A12 = UU(1:L.szsub,(L.szsub+1):end);
    U.A21 = [];
end
