function y = BF_analytic(x,dim)
% compute the analytic counter part of a real signal x
switch dim
    case 1
        yhat = fft(x);
        yhat = yhat*2;
        yhat((1+end/2):end)=0;
        y = ifft(yhat);
    case 2
        yhat = fft2(x);
        yhat = yhat*2;
        yhat((1+end/2):end,:)=0;
        y = ifft2(yhat);
    otherwise
        error('other dimensions are underdevelopment');
end
end
