function X = ifft3(X)
% IFFT3 Three-dimensional inverse discrete Fourier transform.
%    IFFT3(X) returns the three-dimensional inverse Fourier transform of
%    tensor X.
%  
%    See also fft, fft2, fft3, fftn, ifft, ifft2, ifftn.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

X = ifft(X,[],1);
X = ifft(X,[],2);
X = ifft(X,[],3);

end