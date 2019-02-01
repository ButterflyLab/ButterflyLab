function  Factor = HSSBF_update(Afac,Ua,Va,U,V,x,k,rk,tol,mr,lsz,method,pha0)
% HSSBF_update computes the HSSBF of Afac - (Ua*Va).*exp(1i*U*V), where
% Afac is in a HSSBF format.
%
% Input:
% Afac: a structure storing the data-sparse representation of the matrix
%        A(x,k), and some auxiliary variables.
% x and k: denote the row and column indices of the submatrix of A to be
%       inverted and compressed, i.e. invert and compress the submatrix
%       A(x,k).
% Ua, Va, U, and V give a BF matrix (Ua*Va).*exp(1i*U*V)
% r:    the maximum rank of the low-rank approximation
% tol:  the accuracy parameter of the low-rank approximation. If we can
%       achieve a low-rank approximation with a rank smaller than r and an
%       accuracy tol, we will use a smaller rank.
% mr:   the maximum rank of the amplitude and phase recovery
% lsz:  a size parameter; when the matrix A has a size <= lsz, we
%       don't compress
% method: specify which BF to computer HSSBF
%         1: IDBF
%         2: CURBF
%         3: CURBF2
% pha0: phase function in the whole domain
%
% Output:
% Factor: a structure storing the data-sparse representation of the matrix
%        Afac - (Ua*Va).*exp(1i*U*V), and some auxiliary variables.
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

if nargin < 13, pha0 = []; end
if nargin < 12, method = 3; end

M = numel(x); N = numel(k);
Factor = struct('isLeaf',[],'sz',[],'szsub',[],'A11',[],'A12',[],'A21',[],'A22',[]);
if M~=N
    error('matrix is not square');
end
if M>lsz
    Factor.sz = M;
    Factor.szsub = numel(x(1:end/2));
    Factor.isLeaf = 0;
    if isempty(pha0)
        Factor.A11 = HSSBF_update(Afac.A11,Ua(1:end/2,:),Va(:,1:end/2),U(1:end/2,:),V(:,1:end/2),x(1:end/2),k(1:end/2),rk,tol,mr,lsz,method);
        Factor.A22 = HSSBF_update(Afac.A22,Ua((end/2+1):end,:),Va(:,(end/2+1):end),U((end/2+1):end,:),V(:,(end/2+1):end),x((end/2+1):end),k((end/2+1):end),rk,tol,mr,lsz,method);
    else
        pha = @(i,j) pha0(i,j);
        Factor.A11 = HSSBF_update(Afac.A11,Ua(1:end/2,:),Va(:,1:end/2),U(1:end/2,:),V(:,1:end/2),x(1:end/2),k(1:end/2),rk,tol,mr,lsz,method,pha);
        pha = @(i,j) pha0(i+M/2,j+N/2);
        Factor.A22 = HSSBF_update(Afac.A22,Ua((end/2+1):end,:),Va(:,(end/2+1):end),U((end/2+1):end,:),V(:,(end/2+1):end),x((end/2+1):end),k((end/2+1):end),rk,tol,mr,lsz,method,pha);
    end
    
    % fast matvec for A12 - (Ua*Va).*exp(1i*U*V) and its conjugate transpose
    A12afun = @(t,s) (Ua(t,:)*Va(:,s+Factor.szsub)).*exp(1i*U(t,:)*V(:,s+Factor.szsub));
    A12 = @(f) BF_apply(Afac.A12,f);
    A12t = @(f) BF_adj_apply(Afac.A12,f);
    N1 = Afac.szsub; N2 = Afac.sz-Afac.szsub;
    dim = 1;
    % low-rank approximation of the phase function
    if isempty(pha0)
        [AU,AS,AV] = MVBF_Lowrank_Phase(N1,N2,A12,A12t,tol,5*mr,mr,dim,pi,1);
    else
        pha = @(i,j) pha0(i,j+N/2);
        [AU,AS,AV] = MVBF_Lowrank_Phase_initial(N1,N2,A12,A12t,tol,5*mr,mr,dim,pi,pha,1);
    end
    % low-rank approximation of the amplitude function
    [AUa,ASa,AVa] = MVBF_Lowrank_Amp(N1,N2,A12,A12t,tol,5*mr,mr,dim);
    % BF of U12
    AV = AS*AV'; AVa = ASa*AVa';
    funk = @(t,s) (AUa(t,:)*AVa(:,s)).*exp(1i*AU(t,:)*AV(:,s)) - A12afun(t,s);
    tt = 1:N1; tt = tt(:);
    ss = 1:N2; ss = ss(:);
    switch method
        case 1
            Factor.A12 = BF_IDBF(funk,tt,ss,8,rk,tol);
        case 2
            Factor.A12 = CURBF(funk,tt,ss,rk,tol,0);
        case 3
            Factor.A12 = CURBF2(funk,tt,ss,max(rk,8),rk,tol);
    end
    
    % fast matvec for A21 - (Ua*Va).*exp(1i*U*V) and its conjugate transpose
    A21afun = @(t,s) (Ua(t+Factor.szsub,:)*Va(:,s)).*exp(1i*U(t+Factor.szsub,:)*V(:,s));
    A21 = @(f) BF_apply(Afac.A21,f);
    A21t = @(f) BF_adj_apply(Afac.A21,f);
    N1 = Afac.sz-Afac.szsub; N2 = Afac.szsub;
    dim = 1;
    % low-rank approximation of the phase function
    if isempty(pha0)
        [AU,AS,AV] = MVBF_Lowrank_Phase(N1,N2,A21,A21t,tol,5*mr,mr,dim,pi,1);
    else
        pha = @(i,j) pha0(i+M/2,j);
        [AU,AS,AV] = MVBF_Lowrank_Phase_initial(N1,N2,A21,A21t,tol,5*mr,mr,dim,pi,pha,1);
    end
    % low-rank approximation of the amplitude function
    [AUa,ASa,AVa] = MVBF_Lowrank_Amp(N1,N2,A21,A21t,tol,5*mr,mr,dim);
    % BF of U21
    AV = AS*AV'; AVa = ASa*AVa';
    funk = @(t,s) (AUa(t,:)*AVa(:,s)).*exp(1i*AU(t,:)*AV(:,s)) - A21afun(t,s);
    tt = 1:N1; tt = tt(:);
    ss = 1:N2; ss = ss(:);
    switch method
        case 1
            Factor.A21 = BF_IDBF(funk,tt,ss,8,rk,tol);
        case 2
            Factor.A21 = CURBF(funk,tt,ss,rk,tol,0);
        case 3
            Factor.A21 = CURBF2(funk,tt,ss,max(rk,8),rk,tol);
    end
else
    Factor.sz = M;
    Factor.szsub = numel(x(1:end/2));
    Factor.isLeaf = 1;
    K = (Ua*Va).*exp(1i*U*V);
    Factor.A11 = Afac.A11 - K(1:end/2,1:end/2);
    Factor.A22 = Afac.A22 - K((end/2+1):end,(end/2+1):end);
    Factor.A12 = Afac.A12 - K(1:end/2,(end/2+1):end);
    Factor.A21 = Afac.A21 - K((end/2+1):end,1:end/2);
end