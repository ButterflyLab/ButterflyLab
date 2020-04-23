function g = BF_SHT_f2g(f,N)

% compute g(theta,m) from f(psi,theta) using FFT
% f(psi,theta) = \sum_{m=-2N+1}^{2N-1} g(theta,m)exp(i*m*psi)

check = 0;
g = zeros(2*N,4*N-1);
shift = 1/2;
m = -2*N+1:2*N-1;
vec = exp(2*pi*1i*shift/(4*N-1)*m);
if check
    psi = 2*pi*((0:4*N-2)+shift)/(4*N-1);
    mat = exp(1i*psi'*m)'/(4*N-1);
    gext = zeros(2*N,4*N-1);
end
for cnt_theta = 0:2*N-1
    % FFT for the summation in order m
    g(cnt_theta+1,:) = (fftshift(fft(f(:,cnt_theta+1))).')./vec/(4*N-1);
    if check
        gext(cnt_theta+1,:) = (mat*f(:,cnt_theta+1)).';
    end
end
if check
    fprintf('error of stage I: %f\n',norm(gext-g)/norm(gext));
end
end