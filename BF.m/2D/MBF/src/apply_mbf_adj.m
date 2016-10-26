function x = apply_mbf_adj(Factors, y)

if( size(y,2)==0 )
    return;
end

x = zeros(size(y));

for i=1:size(Factors,1)-1
    x(Factors{i,2},:) = apply_bf_adj(Factors{i,1},y);
end
x(Factors{end,2},:) = Factors{end,1}'*y;

end
