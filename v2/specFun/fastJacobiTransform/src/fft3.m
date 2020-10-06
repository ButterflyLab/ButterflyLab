function X = fft3(X)
% FFT3 Three-dimensional discrete Fourier Transform.
%     FFT3(X) returns the three-dimensional Fourier transform of tensor X.
%  
%     See also fft, fft2, fftn, ifft, ifft2, ifft3, ifftn.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

X = fft(X,[],1);
X = fft(X,[],2);
X = fft(X,[],3);

end