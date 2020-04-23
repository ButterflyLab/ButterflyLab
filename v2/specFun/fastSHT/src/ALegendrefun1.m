function mat = ALegendrefun1(x,k,ord,xx,kk,weight)
% x and k are indices
% xx(x) is the root, kk(k) is the degree
% ord is the order
% output a matrix representing the associated Legendre matrix for the given
% order ord, degrees kk(k), and spacial locations xx(x).
% The matrix is generated using the fast algorithms for evaluating the
% Legendre function.

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
        vals = fastALegendre(kk(k(cntd)),ord,xx(x)).*sqrt(weight(x));
        mat(:,cntd) = sgn^ord*vals.';
    else
        mat(:,cntd) = 0;
    end
end
end