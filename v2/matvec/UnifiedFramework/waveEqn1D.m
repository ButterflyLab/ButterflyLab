function [pp,pn,lamp,lamn,Phi] = waveEqn1D(M,N,L,tIdx)
% suppose the time domain is t in [0,1] and M point discretization
% suppose the space domain is x in [0,1] and N point discretization
% solve the equation on a grid of size N by M and interpolate to L by M
% pp,pn,lamp,lamn contain the low-rank factorization of Phi
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

x = (0:1/N:(1-1/N))';
xx = (0:1/L:(1-1/L))';
c = sin(2*pi*x)+2;
[Phip,Phin] = waveEqnPhase1D(M,N,c,1);
Phi = cell(1,numel(tIdx));
pp = cell(1,numel(tIdx));
pn = cell(1,numel(tIdx));
xi = -L/2:1:(L/2-1);
for cnt = 1:numel(tIdx)
    pp{cnt} = BF_fftIntp(Phip(:,tIdx(cnt)),L)+xx;
    pn{cnt} = BF_fftIntp(Phin(:,tIdx(cnt)),L)-xx;
    lamp = abs(xi((1+end/2):end));
    lamn = abs(xi(1:end/2));
    if nargout > 4
        Phi{cnt} = zeros(L,L);
        if 1 % vectorization
            Phi{cnt}(:,1:end/2) = pn{cnt}*lamn;
            Phi{cnt}(:,(1+end/2):end) = pp{cnt}*lamp;
        else
            for cntx = 1:numel(xi)
                if xi(cntx) < 0
                    Phi{cnt}(:,cntx) = abs(xi(cntx))*pn{cnt};
                end
                if xi(cntx) > 0
                    Phi{cnt}(:,cntx) = abs(xi(cntx))*pp{cnt};
                end
            end
        end
    end
end
end