function f = BF_SHT_g2f(g,N)

% reconstructing f(psi,theta) from g(theta,m) using FFT
% f(psi,theta) = \sum_{m=-2N+1}^{2N-1} g(theta,m)exp(i*m*psi)

check = 0;
f = zeros(4*N-1,2*N);
shift = 1/2;
m = -2*N+1:2*N-1;
vec = exp(2*pi*1i*shift/(4*N-1)*m);
if check
    psi = 2*pi*((0:4*N-2)+shift)/(4*N-1);
    mat = exp(1i*psi'*m);
    fext = zeros(4*N-1,2*N);
end
for cnt_theta = 0:2*N-1
    % FFT for the summation in order m
    f(:,cnt_theta+1) = ifft(ifftshift(vec.*g(cnt_theta+1,:)))*(4*N-1).';
    if check
        fext(:,cnt_theta+1) = mat*g(cnt_theta+1,:).';
    end
end
if check
    fprintf('error of stage II: %f\n',norm(fext-f)/norm(fext));
end