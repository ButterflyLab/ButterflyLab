function [A11,A12,A21,A22,flag] = BF_divide1D(A,r,tol,lsz)
% This code divides the butterfly factorization of a matrix A into four
% butterfly factorizations of A11, A12, A21, and A22, where these four
% matrices are the submatrices of A s.t. A=[A11,A12;A21,A22].

if isa(A,'struct')
    % method 1, less coding, more expensive here, but no influence on HSSBF
    M = size(A.U,1); N = size(A.V,2);
    if M >= (2*r)^2 & M >= lsz
        flag = 0;
        kk = 1:N/2; kk = kk(:); kbox = [1,N/2+1];
        xx = 1:M/2; xx = xx(:); xbox = [1,M/2+1];
        
        % first block
        matR = sparse(1:N/2,1:N/2,ones(1,N/2),N,N/2,N/2);
        matL = sparse(1:M/2,1:M/2,ones(1,M/2),M/2,M,M/2);
        Afun = @(f)  matL*BF_apply(A,matR*f);
        Afunt = @(f)  matR'*BF_adj_apply(A,matL'*f);
        A11 = MVBFslow(Afun,Afunt,xx,xbox, kk, kbox, r, tol, 0); % TODO: use a faster method
        
        % second block
        matR = sparse((1:N/2)+N/2,1:N/2,ones(1,N/2),N,N/2,N/2);
        Afun = @(f)  matL*BF_apply(A,matR*f);
        Afunt = @(f)  matR'*BF_adj_apply(A,matL'*f);
        A12 = MVBFslow(Afun,Afunt,xx,xbox, kk, kbox, r, tol, 0); % TODO: use a faster method
        
        % third block
        matR = sparse(1:N/2,1:N/2,ones(1,N/2),N,N/2,N/2);
        matL = sparse(1:M/2,(1:M/2)+M/2,ones(1,M/2),M/2,M,M/2);
        Afun = @(f)  matL*BF_apply(A,matR*f);
        Afunt = @(f)  matR'*BF_adj_apply(A,matL'*f);
        A21 = MVBFslow(Afun,Afunt,xx,xbox, kk, kbox, r, tol, 0); % TODO: use a faster method
        
        % fourth block
        matR = sparse((1:N/2)+N/2,1:N/2,ones(1,N/2),N,N/2,N/2);
        matL = sparse(1:M/2,(1:M/2)+M/2,ones(1,M/2),M/2,M,M/2);
        Afun = @(f)  matL*BF_apply(A,matR*f);
        Afunt = @(f)  matR'*BF_adj_apply(A,matL'*f);
        A22 = MVBFslow(Afun,Afunt,xx,xbox, kk, kbox, r, tol, 0); % TODO: use a faster method
        
        % method 2, more coding to modify BF, very cheap here
        %     A11 = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);
        %     A12 = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);
        %     A21 = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);
        %     A22 = struct('U',[],'GTol',[],'M',[],'HTol',[],'V',[]);
    else
        A = BF_sp2den(A);
    end
end

if ~isa(A,'struct')
    flag = 1;
    A11 = A(1:end/2,1:end/2);
    A12 = A(1:end/2,(1+end/2):end);
    A21 = A((1+end/2):end,1:end/2);
    A22 = A((1+end/2):end,(1+end/2):end);
end
end


