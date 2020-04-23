function mat = ALegendrefun2(x,k,ord,xx,kk,weight)
% x and k are indices
% xx(x) is the root, kk(k) is the degree
% ord is the order
% output a matrix representing the associated Legendre matrix for the given
% order ord, degrees kk(k), and spacial locations xx(x).
% The matrix is generated using the MATLAB-built-in function for evaluating 
% the Legendre function.

if ord < 0
    sgn = -1; 
    ord = -ord; 
else
    sgn = 1;    
end

if isempty(x)
    mat = zeros(0,length(k));
    return;
elseif isempty(k)
    mat = zeros(length(x),0);
    return;
end
mat = zeros(numel(x),numel(k));
for cntd = 1:numel(k)
    if kk(k(cntd)) >= ord
        val = legendre(kk(k(cntd)),cos(xx(x)),'norm');
        mat(:,cntd) = (sgn^ord*val(ord+1,:).').*sqrt(weight(x));
    else
        mat(:,cntd) = 0;
    end
end
end