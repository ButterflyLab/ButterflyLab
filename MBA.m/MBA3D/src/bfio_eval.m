function u = bfio_eval(Nx,Nk,SL,EL,NG,fun,f,mats,dirSL,dirEL)

if EL<2
    error('EL has to be >= 2');
end

grid = bfio_grid(NG);

u = zeros(Nx,Nx,Nx); %?

ML = floor((SL+EL)/2);

nz = 2^EL; %number of zones
zB = Nk/nz;

for z1=0:nz-1
for z2=0:nz-1
for z3=0:nz-1
    if(z1<nz/4 | z2 < nz/4 | z3 < nz/4 | z1 > 3*nz/4-1 | z2 > 3*nz/4-1 | z3 > 3*nz/4-1)
        k1stt = z1*zB;
        k1end = (z1+1)*zB;
        k2stt = z2*zB;
        k2end = (z2+1)*zB;
        k3stt = z3*zB;
        k3end = (z3+1)*zB;

        NOW = [];
        for ell=SL:-1:ML
            nk = 2^ell; %nk: number of boxes in k;
            nx = Nk/nk;  %nx: number of boxes in x; nk*nx = N;
            kB = Nk/nk;  %kB: size of boxes in k;
            xB = 1/nx;  %xB: size of boxes in x; kB*xB = N;

            PRE = NOW;
            NOW = cell(nk,nk,nk);
            for k1=(k1stt/kB):(k1end/kB-1)
            for k2=(k2stt/kB):(k2end/kB-1)
            for k3=(k3stt/kB):(k3end/kB-1)
                NOW{k1+1,k2+1,k3+1} = cell(nx,nx,nx);
                for x1=0:nx-1
                for x2=0:nx-1
                for x3=0:nx-1
                    NOW{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1} = zeros(NG,NG,NG);
                end
                end
                end
            end
            end
            end

            [so1,so2,so3] = ndgrid((0:kB-1)/kB);
            [ko1,ko2,ko3] = ndgrid(grid);

            for k1=(k1stt/kB):(k1end/kB-1)
            for k2=(k2stt/kB):(k2end/kB-1)
            for k3=(k3stt/kB):(k3end/kB-1)
                for x1=0:nx-1
                for x2=0:nx-1
                for x3=0:nx-1
                    xc1 = (x1+1/2)*xB;
                    xc2 = (x2+1/2)*xB;
                    xc3 = (x3+1/2)*xB;
                    %--------------
                    if(ell==SL)
                        s1 = (0:kB-1) + k1*kB - Nk/2;
                        s2 = (0:kB-1) + k2*kB - Nk/2;
                        s3 = (0:kB-1) + k3*kB - Nk/2;
                        %get
                        all = f(s1+Nk/2+1,s2+Nk/2+1,s3+Nk/2+1);
                        %scale
                        trg = [xc1;xc2;xc3];
                        sp1 = (so1+k1)*kB - Nk/2;
                        sp2 = (so2+k2)*kB - Nk/2;
                        sp3 = (so3+k3)*kB - Nk/2;
                        src = [sp1(:)'; sp2(:)'; sp3(:)'];
                        scl = fun(Nk,trg,src);
                        scl = reshape(scl,[kB,kB,kB]);
                        all = all.*scl;
                        %transform
                        tmp = reshape(all, [size(all,1), numel(all)/size(all,1)]);
                        tmp = dirSL*tmp;
                        all = reshape(tmp, [size(tmp,1), size(all,2), size(all,3)]);
                        all = shiftdim(all,1);
                        tmp = reshape(all, [size(all,1), numel(all)/size(all,1)]);
                        tmp = dirSL*tmp;
                        all = reshape(tmp, [size(tmp,1), size(all,2), size(all,3)]);
                        all = shiftdim(all,1);
                        tmp = reshape(all, [size(all,1), numel(all)/size(all,1)]);
                        tmp = dirSL*tmp;
                        all = reshape(tmp, [size(tmp,1), size(all,2), size(all,3)]);
                        all = shiftdim(all,1);
                        %scale
                        trg = [xc1;xc2;xc3];
                        kp1 = (ko1+k1)*kB - Nk/2;
                        kp2 = (ko2+k2)*kB - Nk/2;
                        kp3 = (ko3+k3)*kB - Nk/2;
                        src = [kp1(:)'; kp2(:)'; kp3(:)'];
                        scl = fun(Nk,trg,src);
                        scl = reshape(scl,[NG,NG,NG]);
                        all = all./scl;

                        NOW{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1} = all;
                    else
                        p1 = floor(x1/2);
                        p2 = floor(x2/2);
                        p3 = floor(x3/2);
                        all = zeros(NG,NG,NG);
                        for a=0:1
                        for b=0:1
                        for c=0:1
                            c1 = 2*k1+a;
                            c2 = 2*k2+b;
                            c3 = 2*k3+c;
                            %get
                            cur = PRE{c1+1,c2+1,c3+1}{p1+1,p2+1,p3+1};
                            %scale
                            trg = [xc1;xc2;xc3];
                            cp1 = (ko1+c1)*kB/2 - Nk/2;
                            cp2 = (ko2+c2)*kB/2 - Nk/2;
                            cp3 = (ko3+c3)*kB/2 - Nk/2;
                            src = [cp1(:)'; cp2(:)'; cp3(:)'];
                            scl = fun(Nk,trg,src);
                            scl = reshape(scl,[NG,NG,NG]);
                            cur = cur.*scl;
                            %transform
                            tmp = reshape(cur, [size(cur,1), numel(cur)/size(cur,1)]);
                            tmp = mats{a+1}*tmp;
                            cur = reshape(tmp, [size(tmp,1), size(cur,2), size(cur,3)]);
                            cur = shiftdim(cur,1);
                            tmp = reshape(cur, [size(cur,1), numel(cur)/size(cur,1)]);
                            tmp = mats{b+1}*tmp;
                            cur = reshape(tmp, [size(tmp,1), size(cur,2), size(cur,3)]);
                            cur = shiftdim(cur,1);
                            tmp = reshape(cur, [size(cur,1), numel(cur)/size(cur,1)]);
                            tmp = mats{c+1}*tmp;
                            cur = reshape(tmp, [size(tmp,1), size(cur,2), size(cur,3)]);
                            cur = shiftdim(cur,1);
                            all = all + cur;
                        end
                        end
                        end
                        %scale
                        trg = [xc1;xc2;xc3];
                        kp1 = (ko1+k1)*kB - Nk/2;
                        kp2 = (ko2+k2)*kB - Nk/2;
                        kp3 = (ko3+k3)*kB - Nk/2;
                        src = [kp1(:)'; kp2(:)'; kp3(:)'];
                        scl = fun(Nk,trg,src);
                        scl = reshape(scl,[NG,NG,NG]);
                        all = all./scl;

                        NOW{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1} = all;
                    end
                end
                end
                end
                %x
            end
            end
            end
            %k
            clear PRE;
        end
        clear c1 c2 c3 p1 p2 p3;
        %--------------------------------
        if(1)%switch
            ell = ML;
            nk = 2^ell;
            nx = Nk/nk;
            kB = Nk/nk;
            xB = 1/nx;

            %in unit box
            [ko1,ko2,ko3] = ndgrid(grid);
            [xo1,xo2,xo3] = ndgrid(grid);

            for k1=(k1stt/kB):(k1end/kB-1)
            for k2=(k2stt/kB):(k2end/kB-1)
            for k3=(k3stt/kB):(k3end/kB-1)
                for x1=0:nx-1
                for x2=0:nx-1
                for x3=0:nx-1
                    kp1 = (ko1+k1)*kB - Nk/2;
                    kp2 = (ko2+k2)*kB - Nk/2;
                    kp3 = (ko3+k3)*kB - Nk/2;
                    src = [kp1(:)'; kp2(:)'; kp3(:)'];
                    xp1 = (xo1+x1)*xB;
                    xp2 = (xo2+x2)*xB;
                    xp3 = (xo3+x3)*xB;
                    trg = [xp1(:)'; xp2(:)'; xp3(:)'];
                    all = fun(Nk,trg,src) * NOW{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1}(:);
                    NOW{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1} = reshape(all,[NG,NG,NG]);
                end
                end
                end
            end
            end
            end
        end
        %--------------------------------
        for ell=ML:-1:EL
            nk = 2^ell;
            nx = Nk/nk;
            kB = Nk/nk;
            xB = 1/nx;
            %LLLLLLLLL
            NXT = cell(nk/2,nk/2,nk/2);
            for k1=(k1stt/(2*kB)):(k1end/(2*kB)-1)
            for k2=(k2stt/(2*kB)):(k2end/(2*kB)-1)
            for k3=(k3stt/(2*kB)):(k3end/(2*kB)-1)
                NXT{k1+1,k2+1,k3+1} = cell(2*nx,2*nx,2*nx);
                for x1=0:2*nx-1
                for x2=0:2*nx-1
                for x3=0:2*nx-1
                    NXT{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1} = zeros(NG,NG,NG);
                end
                end
                end
            end
            end
            end

            NT = Nx/nx;
            [to1,to2,to3] = ndgrid((0:NT-1)/NT);
            [xo1,xo2,xo3] = ndgrid(grid);

            for k1=(k1stt/kB):(k1end/kB-1)
            for k2=(k2stt/kB):(k2end/kB-1)
            for k3=(k3stt/kB):(k3end/kB-1)
                kc1 = (k1+1/2)*kB - Nk/2;
                kc2 = (k2+1/2)*kB - Nk/2;
                kc3 = (k3+1/2)*kB - Nk/2;
                for x1=0:nx-1
                for x2=0:nx-1
                for x3=0:nx-1
                    if(ell~=EL)
                        all = NOW{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1};
                        %scale
                        src = [kc1; kc2; kc3];
                        xp1 = (xo1+x1)*xB;
                        xp2 = (xo2+x2)*xB;
                        xp3 = (xo3+x3)*xB;
                        trg = [xp1(:)'; xp2(:)'; xp3(:)'];
                        scl = fun(Nk,trg,src);
                        scl = reshape(scl,[NG,NG,NG]);
                        all = all./scl;
                        %--
                        q1 = floor(k1/2);
                        q2 = floor(k2/2);
                        q3 = floor(k3/2);
                        for a=0:1
                        for b=0:1
                        for c=0:1
                            d1 = 2*x1+a;
                            d2 = 2*x2+b;
                            d3 = 2*x3+c;
                            %transform
                            cur = all;
                            tmp = reshape(cur, [size(cur,1), numel(cur)/size(cur,1)]);
                            tmp = mats{a+1}'*tmp;
                            cur = reshape(tmp, [size(tmp,1), size(cur,2), size(cur,3)]);
                            cur = shiftdim(cur,1);
                            tmp = reshape(cur, [size(cur,1), numel(cur)/size(cur,1)]);
                            tmp = mats{b+1}'*tmp;
                            cur = reshape(tmp, [size(tmp,1), size(cur,2), size(cur,3)]);
                            cur = shiftdim(cur,1);
                            tmp = reshape(cur, [size(cur,1), numel(cur)/size(cur,1)]);
                            tmp = mats{c+1}'*tmp;  cur = reshape(tmp, [size(tmp,1), size(cur,2), size(cur,3)]);
                            cur = shiftdim(cur,1);
                            %scale
                            src = [kc1; kc2; kc3];
                            dp1 = (xo1+d1)*xB/2;
                            dp2 = (xo2+d2)*xB/2;
                            dp3 = (xo3+d3)*xB/2;
                            trg = [dp1(:)'; dp2(:)'; dp3(:)'];
                            scl = fun(Nk,trg,src);
                            scl = reshape(scl,[NG,NG,NG]);
                            cur = cur.*scl;
                            %put
                            NXT{q1+1,q2+1,q3+1}{d1+1,d2+1,d3+1} = NXT{q1+1,q2+1,q3+1}{d1+1,d2+1,d3+1} + cur;
                        end
                        end
                        end
                        %abc
                    else
                        all = NOW{k1+1,k2+1,k3+1}{x1+1,x2+1,x3+1};
                        %scale
                        src = [kc1; kc2; kc3];
                        xp1 = (xo1+x1)*xB;
                        xp2 = (xo2+x2)*xB;
                        xp3 = (xo3+x3)*xB;
                        trg = [xp1(:)'; xp2(:)'; xp3(:)'];
                        scl = fun(Nk,trg,src);
                        scl = reshape(scl,[NG,NG,NG]);
                        all = all./scl;
                        
                        %transform
                        tmp = reshape(all, [size(all,1), numel(all)/size(all,1)]);
                        tmp = dirEL'*tmp;
                        all = reshape(tmp, [size(tmp,1), size(all,2), size(all,3)]);
                        all = shiftdim(all,1);
                        tmp = reshape(all, [size(all,1), numel(all)/size(all,1)]);
                        tmp = dirEL'*tmp;
                        all = reshape(tmp, [size(tmp,1), size(all,2), size(all,3)]);
                        all = shiftdim(all,1);
                        tmp = reshape(all, [size(all,1), numel(all)/size(all,1)]);
                        tmp = dirEL'*tmp;
                        all = reshape(tmp, [size(tmp,1), size(all,2), size(all,3)]);
                        all = shiftdim(all,1);

                        %scale
                        t1 = (0:NT-1) + x1*NT;
                        t2 = (0:NT-1) + x2*NT;
                        t3 = (0:NT-1) + x3*NT;
                        src = [kc1; kc2; kc3];
                        tp1 = (to1+x1)*xB;
                        tp2 = (to2+x2)*xB;
                        tp3 = (to3+x3)*xB;
                        trg = [tp1(:)'; tp2(:)'; tp3(:)'];
                        scl = fun(Nk,trg,src);
                        scl = reshape(scl,[NT,NT,NT]);
                        all = all.*scl;

                        u(t1+1,t2+1,t3+1) = u(t1+1,t2+1,t3+1) + all;
                    end
                end
                end
                end
                %x
            end
            end
            end
            %k
            NOW = NXT;
            clear NXT;
        end
        %----------------
    end
end
end
end
%z


