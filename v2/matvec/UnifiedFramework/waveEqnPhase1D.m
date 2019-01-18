function [Phip,Phin] = waveEqnPhase1D(M,N,c,type)
% suppose the time domain is t in [0,1] and M point discretization
% suppose the space domain is x in [0,1] and N point discretization
% suppose the frequency domain is xi in {-1,1}
% c is a vector of length N
%
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if nargin < 4
    type = 1;
end
if nargin < 3
    c = ones(N,1);
end
if type == -1
    c = -c;
end

dt = 1/M; dx = 1/N;
Phin = zeros(N,M);
Phip = zeros(N,M);
% use the local Lax-Friedrichs Hamiltonian
for cnt = 2:M
    if 0
        if cnt == 2
            vp = Phip(:,cnt-1);
            Phip(:,cnt) = dt*c.*abs((circshift(vp,-1)-circshift(vp,1))/dx/2 -0.5*(circshift(vp,1)+circshift(vp,-1)-2*vp) + 1) + vp;
            vp = Phin(:,cnt-1);
            Phin(:,cnt) = dt*c.*abs((circshift(vp,-1)-circshift(vp,1))/dx/2 -0.5*(circshift(vp,1)+circshift(vp,-1)-2*vp) - 1) + vp;
        else
            vp = Phip(:,cnt-1);
            Phip(:,cnt) = 2*dt*c.*abs((circshift(vp,-1)-circshift(vp,1))/dx/2 -0.5*(circshift(vp,1)+circshift(vp,-1)-2*vp) + 1) + Phip(:,cnt-2);
            vp = Phin(:,cnt-1);
            Phin(:,cnt) = 2*dt*c.*abs((circshift(vp,-1)-circshift(vp,1))/dx/2 -0.5*(circshift(vp,1)+circshift(vp,-1)-2*vp) - 1) + Phip(:,cnt-2);
        end
    else % the third order TVD Runge-Kutta
        vpp = Phip(:,cnt-1);
        vpn = Phin(:,cnt-1);
        u1p = dt*c.*abs((circshift(vpp,-1)-circshift(vpp,1))/dx/2 -0.5*(circshift(vpp,1)+circshift(vpp,-1)-2*vpp) + 1) + vpp;
        u1n = dt*c.*abs((circshift(vpn,-1)-circshift(vpn,1))/dx/2 -0.5*(circshift(vpn,1)+circshift(vpn,-1)-2*vpn) - 1) + vpn;
        u2p = (1/4)*dt*c.*abs((circshift(u1p,-1)-circshift(u1p,1))/dx/2 -0.5*(circshift(u1p,1)+circshift(u1p,-1)-2*u1p) + 1) + (3/4)*vpp + (1/4)*u1p;
        u2n = (1/4)*dt*c.*abs((circshift(u1n,-1)-circshift(u1n,1))/dx/2 -0.5*(circshift(u1n,1)+circshift(u1n,-1)-2*u1n) - 1) + (3/4)*vpn + (1/4)*u1n;
        Phip(:,cnt) = (2/3)*dt*c.*abs((circshift(u2p,-1)-circshift(u2p,1))/dx/2 -0.5*(circshift(u2p,1)+circshift(u2p,-1)-2*u2p) + 1) + (1/3)*vpp + (2/3)*u2p;
        Phin(:,cnt) = (2/3)*dt*c.*abs((circshift(u2n,-1)-circshift(u2n,1))/dx/2 -0.5*(circshift(u2n,1)+circshift(u2n,-1)-2*u2n) - 1) + (1/3)*vpn + (2/3)*u2n;
    end
end