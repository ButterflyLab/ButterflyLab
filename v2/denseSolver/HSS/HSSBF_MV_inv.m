function Factor = HSSBF_MV_inv(Afactor,x,k,r,tol,lsz)
% HSSBF_MV_inv computes the inverse of a HSSBF matrix A stored in a
% data-sparse HSSBF format in the variable Afactor, and returns the inverse
% as a data-sparse HSSBF format in Factor.
%
% Input:
% Afactor: a data-sparse HSSBF format of A.
% x and k: denote the row and column indices of the submatrix of A to be
%       inverted and compressed, i.e. invert and compress the submatrix
%       A(x,k).
% r:    the maximum rank of the low-rank approximation
% tol:  the accuracy parameter of the low-rank approximation. If we can
%       achieve a low-rank approximation with a rank smaller than r and an
%       accuracy tol, we will use a smaller rank.
% lsz:  a size parameter; when the matrix A has a size <= lsz, we
%       directly invert A.
%
% Output:
% Factor: a structure storing the data-sparse representation of the inverse
%       of A(x,k), and some auxiliary variables.
% Factor.isLeaf: 1, if all matrices are dense
%                0, if all matrices are data-sparse
% Factor.sz: the size of A
% Factor.szsub: the size of A11
% Factor.B
% Factor.Btr
% Factor.Bctr
% Factor.A11
% Factor.A22
% Factor.A11inv
% Factor.A22inv
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

%% Step 1: split
M = numel(x); N = numel(k);
Factor = struct('isLeaf',[],'sz',[],'szsub',[],'B',[],'A11',[],'A22',[],'A11inv',[],'A22inv',[]);
if M~=N
    error('matrix is not square');
end
if M>lsz
    M1 = numel(x(1:end/2));
    M2 = M - M1;
    A11_invfac = HSSBF_MV_inv(Afactor.A11,x(1:end/2),k(1:end/2),r,tol,lsz);
    A22_invfac = HSSBF_MV_inv(Afactor.A22,x((end/2+1):end),k((end/2+1):end),r,tol,lsz);
    A12_fac = Afactor.A12;
    A21_fac = Afactor.A21;
    B1_mvfun = @(f) HSSBF_sol(A11_invfac,BF_apply(A12_fac,f));
    B2_mvfun = @(f) HSSBF_sol(A22_invfac,BF_apply(A21_fac,f));
    B1_trmvfun = @(f) BF_trans_apply(A12_fac,HSSBF_trans_sol(A11_invfac,f));
    B2_trmvfun = @(f) BF_trans_apply(A21_fac,HSSBF_trans_sol(A22_invfac,f));
    B1_ctrmvfun = @(f) BF_adj_apply(A12_fac,HSSBF_adj_sol(A11_invfac,f));
    B2_ctrmvfun = @(f) BF_adj_apply(A21_fac,HSSBF_adj_sol(A22_invfac,f));
    % B12 = B1; B21 = B2;
    
    %% Step 2: inverse B22 + I, which is I
    % we do nothing in this step
    
    %% Step 3: compress the Schur complement of B22bar
    % B11ta = -B12*B21;
    % get B11ta in the form of matvec
    B11ta_mvfun = @(f) B1_mvfun(B2_mvfun(-f));
    B11ta_trmvfun = @(f) B1_mvfun(B2_mvfun(-f));
    
    
    B11ta_fun = @(x,k) (Ua(x,:)*Va(:,k)).*exp(1i*U(x,:)*V(:,k));
    
    %% Step 4: inverse I + B11ta
    I_plus_B11ta_fun = @(x,k) (repmat(x(:),[1,M1])==repmat(k(:)',[M1,1])) + B11ta_fun(x,k);
    % compute the inverse of I + B11ta
    I_plus_B11ta_invfac = HSSBF_inv(I_plus_B11ta_fun,1:M1,1:M1,r,tol,lsz);
    % compute B11bar
    B11bar_mvfun = @(f) HSSBF_sol(I_plus_B11ta_invfac,f) - f;
    B11bar_trmvfun = @(f) HSSBF_trans_sol(I_plus_B11ta_invfac,f) - f;
    B11bar_ctrmvfun = @(f) HSSBF_adj_sol(I_plus_B11ta_invfac,f) - f;
    % B11bar = inv(eye(N/2)+B11ta)-eye(N/2);
    
    %% Step 5: compute Bbar
    %Bbarapp_plus_I = [eye(N/2),zeros(N/2);-B21,eye(N/2)]*...
    %    [eye(N/2)+B11bar,zeros(N/2);zeros(N/2),eye(N/2)]*...
    %    [eye(N/2),-B1;zeros(N/2),eye(N/2)];
    Bbar_plus_I_mvfun = @(f) HSSBF_combine(B1_mvfun,B2_mvfun,B11bar_mvfun,M1,f,'n');
    Bbar_plus_I_trmvfun = @(f) HSSBF_combine(B1_trmvfun,B2_trmvfun,B11bar_trmvfun,M1,f,'t');
    Bbar_plus_I_ctrmvfun = @(f) HSSBF_combine(B1_ctrmvfun,B2_ctrmvfun,B11bar_ctrmvfun,M1,f,'c');
    
    % Step 6: inverse A
    %invAapp = (Bbarapp+eye(N))*[inv(A11),zeros(N/2);zeros(N/2),inv(A22)];
    Factor.isLeaf = 0;
    Factor.sz = M;
    Factor.szsub = numel(x(1:end/2));
    Factor.B = Bbar_plus_I_mvfun;
    Factor.Btr = Bbar_plus_I_trmvfun;
    Factor.Bctr = Bbar_plus_I_ctrmvfun;
    Factor.A11inv = A11_invfac;
    Factor.A22inv = A22_invfac;
    Factor.A11 = [];
    Factor.A22 = [];
else
    A11 = Afun(x(1:end/2),k(1:end/2));
    A22 = Afun(x((end/2+1):end),k((end/2+1):end));
    A12 = Afun(x(1:end/2),k((end/2+1):end));
    A21 = Afun(x((end/2+1):end),k(1:end/2));
    B1 = A11\A12;
    B2 = A22\A21;
    % Step 1: split
    B12 = B1; B21 = B2;
    % Step 2: inverse B22 + I, which is I
    % Step 3: compress the Schur complement of B22bar
    B11ta = -B12*B21;
    % Step 4: inverse B11ta
    B11bar = inv(eye(N/2)+B11ta)-eye(N/2);
    % Step 5: compute Bbar
    Bbar_plus_I = [eye(N/2),zeros(N/2);-B21,eye(N/2)]*...
        [eye(N/2)+B11bar,zeros(N/2);zeros(N/2),eye(N/2)]*...
        [eye(N/2),-B1;zeros(N/2),eye(N/2)];
    Factor.isLeaf = 1;
    Factor.sz = M;
    Factor.szsub = numel(x(1:end/2));
    Factor.B = Bbar_plus_I;
    Factor.Btr = Bbar_plus_I.';
    Factor.Bctr = Bbar_plus_I';
    Factor.A11inv = struct('isLeaf',1,'sz',N/2,'szsub',N/4,'B',[],'A11',[],'A22',[],'A11inv',[],'A22inv',[]);
    Factor.A22inv = struct('isLeaf',1,'sz',N/2,'szsub',N/4,'B',[],'A11',[],'A22',[],'A11inv',[],'A22inv',[]);
    Factor.A11 = A11;
    Factor.A22 = A22;
end
end

function y = HSSBF_combine(B1,B2,B11bar,M,p,flavor)
switch flavor
    case 'n'
        y = p;
        y(1:M,:) = y(1:M,:)-B1(p(M+1:end,:));
        y(1:M,:) = y(1:M,:) + B11bar(y(1:M,:));
        y(M+1:end,:) = y(M+1:end,:) - B2(y(1:M,:));
        
    case 't'
        y = p;
        y(1:M,:) = y(1:M,:) - B2(p(M+1:end,:));
        y(1:M,:) = y(1:M,:) + B11bar(y(1:M,:));
        y(M+1:end,:) = y(M+1:end,:)-B1(y(1:M,:));
        
    case 'c'
        y = p;
        y(1:M,:) = y(1:M,:) - B2(p(M+1:end,:));
        y(1:M,:) = y(1:M,:) + B11bar(y(1:M,:));
        y(M+1:end,:) = y(M+1:end,:)-B1(y(1:M,:));
end
end