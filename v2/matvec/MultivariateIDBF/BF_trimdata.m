function d = BF_trimdata(data)
    d = zeros(size(data, 1), 1);
    for l = 1 : size(data, 1)
        tmp = nonzeros(data(l, :));
        d(l) = trimmean(tmp, min(100/length(tmp)+0.1, 49));
    end
end