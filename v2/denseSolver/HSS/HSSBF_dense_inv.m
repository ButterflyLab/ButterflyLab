function invA = HSSBF_dense_inv(A,lsz)
% HSSBF_dense_inv computes the inverse of a HSSBF matrix.
%
% Input:
% A: a HSSBF matrix
% lsz:  a size parameter; when the matrix A has a size <= lsz, we
%       directly invert A.
%
% Output:
% invA
%
%
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Haizhao Yang, 2018, National University of Singapore

%% Step 1: split
[M,N] = size(A);
if M~=N
    error('matrix is not square');
end
if M>lsz
    A11inv = HSSBF_dense_inv(A(1:end/2,1:end/2),lsz);
    A22inv = HSSBF_dense_inv(A((end/2+1):end,(end/2+1):end),lsz);
    A12 = A(1:end/2,(end/2+1):end);
    A21 = A((end/2+1):end,1:end/2);
    B1 = A11inv*A12; % B1 is available in a form of fast matvec
    B2 = A22inv*A21; % B2 is available in a form of fast matvec
    
    %% Step 2: inverse B22 + I, which is I
    % we do nothing in this step
    
    %% Step 3: compress the Schur complement of B22bar
    B11ta = -B1*B2; % B11ta is available in a form of fast matvec
    
    %% Step 4: inverse I + B11ta
    B11ta_plus_I = B11ta+eye(N/2); % B11ta_plus_I is available in a form of fast matvec
    B11bar = HSSBF_dense_inv(B11ta_plus_I,lsz)-eye(N/2); % B11bar is available in a form of fast matvec
    
    %% Step 5: compute Bbar
    Bbarapp_plus_I = [eye(N/2),zeros(N/2);-B2,eye(N/2)]*...
       [eye(N/2)+B11bar,zeros(N/2);zeros(N/2),eye(N/2)]*...
       [eye(N/2),-B1;zeros(N/2),eye(N/2)];
    % Bbarapp_plus_I is available in a form of fast matvec
    % Step 6: inverse A
    invA = Bbarapp_plus_I*[A11inv,zeros(N/2);zeros(N/2),A22inv];
    % invA is available in a form of fast matvec
else
    A11 = A(1:end/2,1:end/2);
    A22 = A((end/2+1):end,(end/2+1):end);
    A12 = A(1:end/2,(end/2+1):end);
    A21 = A((end/2+1):end,1:end/2);
    B1 = A11\A12;
    B2 = A22\A21;
    % Step 1: split
    % Step 2: inverse B22 + I, which is I
    % Step 3: compress the Schur complement of B22bar
    B11ta = -B1*B2;
    % Step 4: inverse B11ta
    B11bar = inv(eye(N/2)+B11ta)-eye(N/2);
    % Step 5: compute Bbar
    Bbarapp_plus_I = [eye(N/2),zeros(N/2);-B2,eye(N/2)]*...
        [eye(N/2)+B11bar,zeros(N/2);zeros(N/2),eye(N/2)]*...
        [eye(N/2),-B1;zeros(N/2),eye(N/2)];
    % Step 6: inverse A
    invA = Bbarapp_plus_I*[inv(A11),zeros(N/2);zeros(N/2),inv(A22)];
end
end
