function u = BF_fftIntp(v,N)
v = v(:);
n = numel(v);
vh = fftshift(fft(v));
vhext = zeros(N,1);
vhext(((N-n)/2+1):(end-(N-n)/2)) = vh;
u = real(ifft(ifftshift(vhext))*sqrt(N/n))*2;
end