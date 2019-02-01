function Factor = HSSBF_RS_inv(Afun,x,k,r,tol,lsz)
% HSSBF_RS_inv computes the inverse of a HSSBF matrix stored in a
% function handle Afun and returns the inverse as a data-sparse HSSBF
% format in Factor.
%
% Input:
% Afun: a function handle for evaluating an arbitrary entry of the HSSBF
%       matrix A in O(1) operations, i.e. A(i,j) = Afun(i,j).
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
    A11_invfac = HSSBF_RS_inv(Afun,x(1:end/2),k(1:end/2),r,tol,lsz);
    A22_invfac = HSSBF_RS_inv(Afun,x((end/2+1):end),k((end/2+1):end),r,tol,lsz);
    A12_fac = CURBF(Afun,x(1:end/2),k((end/2+1):end),r,tol);
    A21_fac = CURBF(Afun,x((end/2+1):end),k(1:end/2),r,tol);
    B1_mvfun = @(f) HSSBF_sol(A11_invfac,BF_apply(A12_fac,f));
    B2_mvfun = @(f) HSSBF_sol(A22_invfac,BF_apply(A21_fac,f));
    B1_trmvfun = @(f) BF_trans_apply(A12_fac,HSSBF_trans_sol(A11_invfac,f));
    B2_trmvfun = @(f) BF_trans_apply(A21_fac,HSSBF_trans_sol(A22_invfac,f));
    B1_ctrmvfun = @(f) BF_adj_apply(A12_fac,HSSBF_adj_sol(A11_invfac,f));
    B2_ctrmvfun = @(f) BF_adj_apply(A21_fac,HSSBF_adj_sol(A22_invfac,f));
    % B12 = B1; B21 = B2;
    if 0 % for debug only, exact B1 and B2
        A11 = Afun(x(1:end/2),k(1:end/2));
        A22 = Afun(x((end/2+1):end),k((end/2+1):end));
        A12 = Afun(x(1:end/2),k((end/2+1):end));
        A21 = Afun(x((end/2+1):end),k(1:end/2));
        B1 = A11\A12;
        B2 = A22\A21;
        B1_mvfun = @(f) B1*f;
        B2_mvfun = @(f) B2*f;
        B1_trmvfun = @(f) B1.'*f;
        B2_trmvfun = @(f) B2.'*f;
        B1_ctrmvfun = @(f) B1'*f;
        B2_ctrmvfun = @(f) B2'*f;
    end
    
    %% Step 2: inverse B22 + I, which is I
    % we do nothing in this step
    
    %% Step 3: compress the Schur complement of B22bar
    % B11ta = -B1*B2;
    
    if 0
        % compute the low-rank information of B11ta
        if 0 % no mode decomposition
            % get B11ta in the form of matvec
            B11ta_mvfunr = @(f) real(B1_mvfun(B2_mvfun(-f)));
            B11ta_mvfuni = @(f) imag(B1_mvfun(B2_mvfun(-f)));
            B11ta_trmvfunr = @(f) real(B1_mvfun(B2_mvfun(-f)));
            B11ta_trmvfuni = @(f) imag(B1_mvfun(B2_mvfun(-f)));
            % low-rank approximation of the amplitude function
            [Uar,Sar,Var] = MVBF_Lowrank_Amp(M/2,N/2,B11ta_mvfunr,B11ta_trmvfunr,tol,r+2,r,1,1);
            [Uai,Sai,Vai] = MVBF_Lowrank_Amp(M/2,N/2,B11ta_mvfuni,B11ta_trmvfuni,tol,r+2,r,1,1);
            % low-rank approximation of the phase function
            [Ur,Sr,Vr] = MVBF_Lowrank_Phase(M/2,N/2,B11ta_mvfunr,B11ta_trmvfunr,tol,r+2,r,1,pi*2,0,1);
            [Ui,Si,Vi] = MVBF_Lowrank_Phase(M/2,N/2,B11ta_mvfuni,B11ta_trmvfuni,tol,r+2,r,1,pi*2,0,1);
            
            % reconstruct the B11ta matrix
            Vr = Sr*Vr'; Vi = Si*Vi'; Var = Sar*Var'; Vai = Sai*Vai';
            B11ta_fun = @(x,k) complex(  (Uar(x,:)*Var(:,k)).*(cos(Ur(x,:)*Vr(:,k))),  (Uai(x,:)*Vai(:,k)).*( cos(Ui(x,:)*Vi(:,k)) ) );
        else % use mode decomposition
            % get B11ta in the form of matvec
            B11ta_mvfun = @(f) B1_mvfun(B2_mvfun(-f));
            B11ta_trmvfun = @(f) B1_mvfun(B2_mvfun(-f));
            opt.numCom = 3;
            opt.VMDonly = 1;
            [U,S,V,Ua,Sa,Va] = MVBF_Lowrank_MD(M/2,N/2,B11ta_mvfun,B11ta_trmvfun,tol,r+2,r,1,pi*2,0,1,opt);
            
            % reconstruct the B11ta matrix
            B11ta_funCell = cell(1,opt.numCom);
            for cntc = 1:opt.numCom
            V{cntc} = S{cntc}*V{cntc}'; Va{cntc} = Sa{cntc}*Va{cntc}'; 
            B11ta_funCell{cntc} = @(x,k) (Ua{cntc}(x,:)*Va{cntc}(:,k)).*(exp(1i*U{cntc}(x,:)*V{cntc}(:,k)));
            end
            B11ta_fun = @(x,k) B11ta_funCell{1}(x,k)+B11ta_funCell{2}(x,k)+B11ta_funCell{3}(x,k);
        end
        %         B11taapp = B11ta_fun(eye(N/2));
        %         B11ta = B11ta_mvfun(eye(N/2));
        %         if size(B11ta,1)>256
        %             save('data.mat','B11ta','B11taapp');
        %         end
        %         B11ta_fun = @(x,k) B11ta(x,k);
    else % for the purpose of debug only
        B11ta_mvfun = @(f) B1_mvfun(B2_mvfun(-f));
        B11ta = B11ta_mvfun(eye(N/2));
        
        if 0
            % test the influence of VMD
            diff = BF_VMD(B11ta);
            %            display('residual error');
            %            norm(diff)/norm(B11ta)
            
            B11ta = B11ta - diff;
        end
        %         if size(B11ta,1)>256
        %             save('data.mat','B11ta');
        %         end
        B11ta_fun = @(x,k) B11ta(x,k);
    end
    
    %% Step 4: inverse I + B11ta
    I_plus_B11ta_fun = @(x,k) (repmat(x(:),[1,numel(k)])==repmat(k(:)',[numel(x),1])) + B11ta_fun(x,k);
    % compute the inverse of I + B11ta
    I_plus_B11ta_invfac = HSSBF_RS_inv(I_plus_B11ta_fun,1:M1,1:M1,r,tol,lsz);
    % compute B11bar
    B11bar_mvfun = @(f) HSSBF_sol(I_plus_B11ta_invfac,f) - f;
    B11bar_trmvfun = @(f) HSSBF_trans_sol(I_plus_B11ta_invfac,f) - f;
    B11bar_ctrmvfun = @(f) HSSBF_adj_sol(I_plus_B11ta_invfac,f) - f;
    %             B11bar = inv(eye(N/2)+B11ta)-eye(N/2);
    
    %% Step 5: compute Bbar
    %     Bbarapp_plus_I = [eye(N/2),zeros(N/2);-B2,eye(N/2)]*...
    %        [eye(N/2)+B11bar,zeros(N/2);zeros(N/2),eye(N/2)]*...
    %        [eye(N/2),-B1;zeros(N/2),eye(N/2)];
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
    Bbar_plus_I = [eye(N/2),zeros(N/2);-B2,eye(N/2)]*...
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