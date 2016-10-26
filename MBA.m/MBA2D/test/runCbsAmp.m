function runCbsAmp(N, func_name, EPS, fig)

addpath('../src/');

EL = 3;
SL = log2(N)-EL;
stoplev = 5;

switch func_name
    case 'funF'
        fun = @(x,k)funF(N,x,k);
    case 'fun0'
        fun = @(x,k)fun0(N,x,k);
        aun = @(x,k)aun0(N,x,k);
    case 'fun1'
        fun = @(x,k)fun1(N,x,k);
        aun = @(x,k)aun1(N,x,k);
end

[mats,dir,dirlev] = bfio_prep(EL,EPS,N,stoplev);
[arp,~] = bfio_setup(N, aun);

if(1)
    f = randn(N,N) + sqrt(-1)*randn(N,N);
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'w');
    string = {'CpxNumMat'};
    serialize(fid, f, string);
end
if(1)
    binstr = sprintf('f_%d.bin', N);
    fid = fopen(binstr,'r');
    string = {'CpxNumMat'};
    f = deserialize(fid, string);
end


maskcase = 4;
switch maskcase
    case 1
        mask =ones(N,N);
        mask(end/4+1:3*end/4,end/4+1:3*end/4) = 0;
    case 2
        mask =ones(N,N);
        mask(3*end/8+1:5*end/8,3*end/8+1:5*end/8) = 0;
    case 3
        mask =ones(N,N);
        mask(end/4+1:3*end/4,end/4+1:3*end/4) = 0;
        mask(7*end/16+1:9*end/16,7*end/16+1:9*end/16) = 1;
    case 4
        mask = ones(N,N);
end
f = f.*mask;

t0 = cputime;
u = bfioChebyshevAmp(N,N,SL,EL,EPS,fun,aun,arp,f,mats,dir,dirlev,stoplev,1);
te = cputime-t0;

NC = 256;
t0 = cputime;
relerr = bfio_check_amp(N,fun,aun,f,u,NC);
tc = (cputime-t0)*N*N/NC;
rt = tc/te;

fprintf(fig,'---------------------------------\n');
fprintf(fig,'N          : %4d\n', N);
fprintf(fig,'EPS        : %4d\n', EPS);
fprintf(fig,'stoplev    : %4d\n', stoplev);
fprintf(fig,'relerr     : %.3e\n', relerr);
fprintf(fig,'eval time  : %.3e\n',te);
fprintf(fig,'check time : %.3e\n',tc);
fprintf(fig,'ratio      : %.3e\n', rt);
fprintf(fig,'---------------------------------\n\n');

end
