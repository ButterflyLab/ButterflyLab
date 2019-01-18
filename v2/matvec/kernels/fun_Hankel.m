function res = fun_Hankel(N,x,k)
% Bessel function in MATLAB

x = N + 2*pi/3*(x-1);
k = k-1;
res = besselh(repmat(k',length(x),1),repmat(x,1,length(k)));

end
