function fun = directjac3d(nts,r,s,t,wghtr,wghts,wghtt,n,da,db,c,flag)
if nargin < 10
    flag = [];
end

if nts < 2^12
   it = 10;
else
   it = 28;
end
valr = jacrecur(nts,r,wghtr,it-1,da,db);
vals = jacrecur(nts,s,wghts,it-1,da,db);
valt = jacrecur(nts,t,wghtt,it-1,da,db);
m = length(n);

nu = [it:nts-1]';
nt = zeros(nts,1);

if  flag ~= 'FULL'
    fun = @(c)dir3d1(c)

else
    Jr = interpjac1(nt,r,nu,da,db,-1);
    Jr = real(Jr.*exp(1i*r*nu'));
    Jr = [valr Jr].';
    Js = interpjac1(nt,s,nu,da,db,-1);
    Js = real(Js.*exp(1i*s*nu'));
    Js = [vals Js].';
    Jt = interpjac1(nt,t,nu,da,db,-1);
    Jt = real(Jt.*exp(1i*t*nu'));
    Jt = [valt Jt];
    fun = @(c)dir3d2(c);

end
    function w = dir3d1(c)
        c = c(:);
        w = zeros(m,1);
        for i = 1:m

            kk = floor((n(i)-0.5)/nts/nts);
            jj = floor((n(i)-kk*nts*nts-0.5)/nts);
            ii = n(i)-kk*nts*nts-nts*jj;
            x = interpjac1(nt,r(kk+1),nu,da,db,-1);
            %x1 = exp(1i*2*pi/nts*floor(s(jj+1)*nts/2/pi)*nu');
            x1 = exp(1i*r(kk+1)*nu');
            x = real(x.*x1);
            y = interpjac1(nt,s(jj+1),nu,da,db,-1);
            %z1 = exp(1i*2*pi/nts*floor(t(ii)*nts/2/pi)*nu');
            y1 = exp(1i*s(jj+1)*nu');
            y = real(y.*y1);
            z = interpjac1(nt,t(ii),nu,da,db,-1);
            %z1 = exp(1i*2*pi/nts*floor(t(ii)*nts/2/pi)*nu');
            z1 = exp(1i*t(ii)*nu');
            z = real(z.*z1);
            x = [valr(kk+1,:) x];
            y = [vals(jj+1,:) y];
            z = [valt(ii,:) z];
            w(i) = real(kron(kron(x,y),z)*c);
        end
    end

    function y = dir2d(c)
        y = Jt*c*Js;
        y = y(:);
    end
    function y = dir3d2(c)
        
        y = zeros(nts^2,nts);
        for i = 1:nts
            y(:,i) = dir2d(c(:,:,i));
        end
        y = y*Jr;
        y = y(:);
    end
end
