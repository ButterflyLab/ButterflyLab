close all;
addpath('../src/');

fi = 1;

N = 256;
EPS = 5;
EL = 3;
SL = log2(N)-EL;
stoplev = 5;
if(1)
    switch fi
        case 0
            fun = @(x,k)fun0(N,x,k);
        case 1
            fun = @(x,k)fun1(N,x,k);
    end
    [mats,dir,dirlev] = bfio_prep(EL,EPS,N,stoplev);
    
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
    u = bfioChebyshev(N,N,SL,EL,EPS,fun,f,mats,dir,dirlev,stoplev,1);
    te = cputime-t0;
    
    NC = 256;
    t0 = cputime;
    relerr = bfio_check(N,fun,f,u,NC);
    tc = (cputime-t0)*N*N/NC;
    rt = tc/te;
    
    fprintf(1,'N %d\n', N);
    fprintf(1,'EPS %d\n', EPS);
    fprintf(1,'EL %d\n', EL);
    fprintf(1,'stoplev %d\n', stoplev);
    fprintf(1,'relerr %d\n', relerr);
    fprintf(1,'eval time %d\n',te);
    fprintf(1,'check time %d\n',tc);
    fprintf(1,'ratio %d\n', rt);
end



