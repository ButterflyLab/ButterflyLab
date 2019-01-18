function y = fun_Bessel(nu,r,w)
% This function evaluate the Bessel function of order nu at the point w*r
%
% Copyright 2017 by Yingzhou Li and Haizhao Yang

y = besselj(nu,r*w');
end