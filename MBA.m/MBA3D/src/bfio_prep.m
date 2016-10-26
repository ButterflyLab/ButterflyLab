function [mats,dir,dirlev] = bfio_prep(EL,EPS,N,stoplev)
    if N>2^stoplev
        gs = bfio_grid(EPS);

        mats = cell(2,1);
        ts = gs/2;
        mats{1} = bfio_prep_aux(gs,ts);
        ts = gs/2+1/2;
        mats{2} = bfio_prep_aux(gs,ts);

        NT = 2^EL;
        ts = [0:NT-1]/NT;
        dir = bfio_prep_aux(gs,ts);

        numlev = log2(N)-stoplev;
        dirlev = cell(1,numlev);
        for lev = 1:numlev
            NT = 2^(EL+lev-1);
            ts = [0:NT-1]/NT;
            dirlev{lev} = bfio_prep_aux(gs,ts);
        end
    else
        dir = [];
        mats = [];
        dirlev = [];
    end
end