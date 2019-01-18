function y = apply_maskbf(maskFactor, x)

y = zeros(size(x));
for it = 1:length(maskFactor)
    if maskFactor{it}.isfactor
        factor = maskFactor{it}.factor;
        xidx = maskFactor{it}.xidx;
        kidx = maskFactor{it}.kidx;
        y(xidx,:) = y(xidx,:) + BF_apply(factor, x(kidx,:));
    else
        factor = maskFactor{it}.factor;
        xidx = maskFactor{it}.xidx;
        kidx = maskFactor{it}.kidx;
        y(xidx,:) = y(xidx,:) + factor*x(kidx,:);
    end
end

end