
% Reference: A Hierarchical Butterfly LU Preconditioner for Two-Dimensional
% Electromagnetic Scattering Problems Involving Open Surfaces, preprint,
% 2019
% Copyright reserved.
% By Yang Liu and Haizhao Yang, 2018

origin = [200,60];
ppw = 15;

switch num
    case 1
        %% cycle
        a = 1.0;
        b = 2.0;
        st = 0;
        ed = 2*pi;%1.95 not OK
        t = linspace(st,ed,N);
        dt = (ed-st)/N;
        rho1 = [a*cos(t-dt/2) ; b*sin(t-dt/2)];
        rho2 = [a*cos(t+dt/2) ; b*sin(t+dt/2)];
        dl = sqrt(sum((rho2-rho1).^2));
        rho = [a*cos(t) ; b*sin(t)];
    case 2
        %% spiral lines
        a = 1.0;
        b = 1.0;
        st = 0;
        ed = 2*pi;
        t = linspace(st,ed,N);
        dt = (ed-st)/N;
        r_st = 1.0;
        r_ed = 2.0;
        rad = linspace(r_st,r_ed,N);
        drad = (r_ed-1)/N;
        rho1 = [a*cos(t-dt/2).*(rad-drad/2) ; b*sin(t-dt/2).*(rad-drad/2)];
        rho2 = [a*cos(t+dt/2).*(rad+drad/2) ; b*sin(t+dt/2).*(rad+drad/2)];
        dl = sqrt(sum((rho2-rho1).^2));
        rho = [a*cos(t).*rad ; b*sin(t).*rad];
    case 3
        %% half cycle
        a = 1.0;
        b = 2.0;
        st = 0;
        ed = pi;
        t = linspace(st,ed,N);
        dt = (ed-st)/N;
        rho1 = [a*cos(t-dt/2) ; b*sin(t-dt/2)];
        rho2 = [a*cos(t+dt/2) ; b*sin(t+dt/2)];
        dl = sqrt(sum((rho2-rho1).^2));
        rho = [a*cos(t) ; b*sin(t)];
    case 4
        %% spiral lines
        a = 1.0;
        b = 1.0;
        st = 0;
        ed = 0.6*pi;
        t = linspace(st,ed,N);
        dt = (ed-st)/N;
        r_st = 1.0;
        r_ed = 2.0;
        rad = linspace(r_st,r_ed,N);
        drad = (r_ed-1)/N;
        rho1 = [a*cos(t-dt/2).*(rad-drad/2) ; b*sin(t-dt/2).*(rad-drad/2)];
        rho2 = [a*cos(t+dt/2).*(rad+drad/2) ; b*sin(t+dt/2).*(rad+drad/2)];
        dl = sqrt(sum((rho2-rho1).^2));
        rho = [a*cos(t).*rad ; b*sin(t).*rad];
        
    case 5
        %% parallel strip
        N = N - 1;
        rho = zeros(2,N);
        rho(1,1:end/2)=-1;
        rho(2,1:end/2)=linspace(-1,1,N/2);
        rho(1,end/2+1:end)=1;
        rho(2,end/2+1:end)=linspace(-1,1,N/2);
        dl = ones(1,N)*4/N;
        
    case 6
        %% rectangular cup
        N = N - 1;
        rho = [ones(1,N/4)*(-2), linspace(-2,2,N/2),ones(1,N/4)*2];
        rho = [rho;linspace(-1,1,N/4), zeros(1,N/2), linspace(-1,1,N/4)];
        dl = ones(1,N)*8/N;
        
    case 7
        %% corrugated line
        a = 1.0;
        b = 1;
        st = 0;
        ed = pi;
        Nseg=3;
        N = N + 2;
        t = linspace(st,ed,N/Nseg);
        dt = (ed-st)/(N/Nseg);
        rho1 = [a*cos(t-dt/2) ; b*sin(t-dt/2)];
        rho2 = [a*cos(t+dt/2) ; b*sin(t+dt/2)];
        dl0 = sqrt(sum((rho2-rho1).^2));
        rho0 = [a*cos(t) ; b*sin(t)];
        rho = [];
        dl=[];
        s=-1;
        for ii=1:Nseg
            rho0(2,:)=rho0(2,:)*s;
            rho0(1,:)=rho0(1,:)+2*a;
            rho = [rho,rho0];
            dl = [dl, dl0];
        end
    case 8
        %% kite
        st = 0;
        ed = 2*pi;
        t = linspace(st,ed,N);
        dt = (ed-st)/N;
        a = 2.0+ sin(3*t);
        b = 2.0+ sin(3*t);
        rho1 = [a.*cos(t-dt/2) ; b.*sin(t-dt/2)];
        rho2 = [a.*cos(t+dt/2) ; b.*sin(t+dt/2)];
        dl = sqrt(sum((rho2-rho1).^2));
        rho = [a.*cos(t) ; b.*sin(t)];
    case 9
        %% bean
        st = 0;
        ed = 2*pi;
        t = linspace(st,ed,N);
        dt = (ed-st)/N;
        a = 2.0+ sin(2*t);
        b = 2.0+ sin(2*t);
        rho1 = [a.*cos(t-dt/2) ; b.*sin(t-dt/2)];
        rho2 = [a.*cos(t+dt/2) ; b.*sin(t+dt/2)];
        dl = sqrt(sum((rho2-rho1).^2));
        rho = [a.*cos(t) ; b.*sin(t)];
    case 10
        %% partial kite
        st = 0;
        ed = 1.8*pi;
        t = linspace(st,ed,N);
        dt = (ed-st)/N;
        a = 2.0+ sin(3*t);
        b = 2.0+ sin(3*t);
        rho1 = [a.*cos(t-dt/2) ; b.*sin(t-dt/2)];
        rho2 = [a.*cos(t+dt/2) ; b.*sin(t+dt/2)];
        dl = sqrt(sum((rho2-rho1).^2));
        rho = [a.*cos(t) ; b.*sin(t)];
end
figure;plot(rho(1,:),rho(2,:),'.');
p = sum(dl);
lambda = p/N*ppw;
k = 2*pi/lambda;
mu0 = 4*pi*1.0e-7;
eps0 = 8.854187*1.0e-12;
c = 1/sqrt(mu0*eps0);
omega = 2*pi*c/lambda;
gamma = 1.781072418;

sc = 1;
Z0fun = @(i,j) HSSBF_Zfun(i,j,num,N,omega,mu0,dl,k,rho,gamma,sc);
Zampfun = @(i,j) HSSBF_Zampfun(i,j,num,N,omega,mu0,dl,k,rho,gamma,sc);
Zpha = @(i,j) HSSBF_Zphafun(i,j,num,N,omega,mu0,dl,k,rho,gamma,sc);
switch num
    case 5
        N = N;
    case 6
        N = N;
    case 7
        N = N-3;
    otherwise
        N = N - 1;
end
nZ = max(abs(Z0fun(1:N,1)));
Z = @(i,j) Z0fun(i,j)/nZ;
Zamp = @(i,j) Zampfun(i,j)/nZ;
% end
f = randn(N,1) + sqrt(-1)*randn(N,1);

