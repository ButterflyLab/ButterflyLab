function Factor = HSSTP_inv(HTP,N,tol,lsz)
% HSSTP_inv computes the inverse of a HSSTP matrix stored in a hierarchical
% Toeplitz structure HTP and returns the inverse as a data-sparse HSSTP
% format in Factor.
%
% Input:
% HTP:  a structure containing four parts, A11, A12r, A12c, A21r, A21c, and
%       A22, where A11 and A22 are two sub HTP, A12r and A12c are the first
%       row and column of A12 block, and similarly for A21r and A21c.
% N:    size of A.
%       first column of the degree of freedom in the off-diagonal blocks.
% tol:  the accuracy parameter of the low-rank approximation. If we can
%       achieve a low-rank approximation with a rank smaller than r and an
%       accuracy tol, we will use a smaller rank.
% lsz:  a size parameter; when the matrix A has a size <= lsz, we
%       directly invert A.
%
% Output:
% Factor: a structure storing the data-sparse representation of the inverse
%       of A, and some auxiliary variables.
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

Factor = struct('isLeaf',[],'sz',[],'szsub',[],'B',[],'A11',[],'A22',[],'A11inv',[],'A22inv',[]);
if N>lsz
    A11_invfac = HSSTP_inv(HTP,N/2,tol,lsz);
    A22_invfac = HSSTP_inv(HTP,N/2,tol,lsz);
    B1_mvfun = @(f) HSSBF_sol(A11_invfac,TP_apply1(A12r,A12c,f));
    B2_mvfun = @(f) HSSBF_sol(A22_invfac,TP_apply1(A21r,A21c,f));
    B1_trmvfun = @(f) TP_trans_apply1(A12r,A12c,HSSBF_trans_sol(A11_invfac,f));
    B2_trmvfun = @(f) TP_trans_apply(A12r,A12c,HSSBF_trans_sol(A22_invfac,f));
    B1_ctrmvfun = @(f) TP_adj_apply1(A12r,A12c,HSSBF_adj_sol(A11_invfac,f));
    B2_ctrmvfun = @(f) TP_adj_apply1(A12r,A12c,HSSBF_adj_sol(A22_invfac,f));
    
    % get B11ta in the form of matvec, B11ta is a Toeplitz matrix
    B11ta_mvfun = @(f) B1_mvfun(B2_mvfun(-f));
    B11ta_trmvfun = @(f) B1_mvfun(B2_mvfun(-f));
     
    % get the first row and column of B11ta
    vec = zeros(N/2,1); vec(1) = 1;
    B11tar = B11ta_trmvfun(vec);
    B11tac = B11ta_mvfun(vec);
    
    % inverse I + B11ta, which is a Toeplitz matrix
    B11tar(1) = B11tar(1)+1; B11tac(1) = B11tac(1)+1;
    B11tar_HTP = 
    I_plus_B11ta_invfac = HSSTP_inv(I_plus_B11ta_fun,1:M1,1:M1,r,tol,lsz);
    % compute B11bar
    B11bar_mvfun = @(f) HSSTP_sol(I_plus_B11ta_invfac,f) - f;
    B11bar_trmvfun = @(f) HSSBF_trans_sol(I_plus_B11ta_invfac,f) - f;
    B11bar_ctrmvfun = @(f) HSSTP_adj_sol(I_plus_B11ta_invfac,f) - f;
    %             B11bar = inv(eye(N/2)+B11ta)-eye(N/2);
    
    %% Step 5: compute Bbar
    %     Bbarapp_plus_I = [eye(N/2),zeros(N/2);-B2,eye(N/2)]*...
    %        [eye(N/2)+B11bar,zeros(N/2);zeros(N/2),eye(N/2)]*...
    %        [eye(N/2),-B1;zeros(N/2),eye(N/2)];
    Bbar_plus_I_mvfun = @(f) HSSTP_combine(B1_mvfun,B2_mvfun,B11bar_mvfun,M1,f,'n');
    Bbar_plus_I_trmvfun = @(f) HSSTP_combine(B1_trmvfun,B2_trmvfun,B11bar_trmvfun,M1,f,'t');
    Bbar_plus_I_ctrmvfun = @(f) HSSTP_combine(B1_ctrmvfun,B2_ctrmvfun,B11bar_ctrmvfun,M1,f,'c');
    
    % Step 6: inverse A
    %invAapp = (Bbarapp+eye(N))*[inv(A11),zeros(N/2);zeros(N/2),inv(A22)];
    Factor.isLeaf = 0;
    Factor.sz = N;
    Factor.szsub = numel(x(1:end/2));
    Factor.B = Bbar_plus_I_mvfun;
    Factor.Btr = Bbar_plus_I_trmvfun;
    Factor.Bctr = Bbar_plus_I_ctrmvfun;
    Factor.A11inv = A11_invfac;
    Factor.A22inv = A22_invfac;
    Factor.A11 = [];
    Factor.A22 = [];
else
    A11 = HTP(x(1:end/2),k(1:end/2));
    A22 = HTP(x((end/2+1):end),k((end/2+1):end));
    A12 = HTP(x(1:end/2),k((end/2+1):end));
    A21 = HTP(x((end/2+1):end),k(1:end/2));
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
    Bbar_plus_I = [eye(N/2),zeros(N/2);-B2,eye(N/2)]*...
        [eye(N/2)+B11bar,zeros(N/2);zeros(N/2),eye(N/2)]*...
        [eye(N/2),-B1;zeros(N/2),eye(N/2)];
    Factor.isLeaf = 1;
    Factor.sz = N;
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

function y = HSSTP_combine(B1,B2,B11bar,N,p,flavor)
switch flavor
    case 'n'
        y = p;
        y(1:N,:) = y(1:N,:)-B1(p(N+1:end,:));
        y(1:N,:) = y(1:N,:) + B11bar(y(1:N,:));
        y(N+1:end,:) = y(N+1:end,:) - B2(y(1:N,:));
        
    case 't'
        y = p;
        y(1:N,:) = y(1:N,:) - B2(p(N+1:end,:));
        y(1:N,:) = y(1:N,:) + B11bar(y(1:N,:));
        y(N+1:end,:) = y(N+1:end,:)-B1(y(1:N,:));
        
    case 'c'
        y = p;
        y(1:N,:) = y(1:N,:) - B2(p(N+1:end,:));
        y(1:N,:) = y(1:N,:) + B11bar(y(1:N,:));
        y(N+1:end,:) = y(N+1:end,:)-B1(y(1:N,:));
end
end