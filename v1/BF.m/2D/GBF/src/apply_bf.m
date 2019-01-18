function y = apply_bf(Factor, y)

if( size(y,2)==0 )
    return;
end
y = Factor.V'*y;

for i=length(Factor.BTol):-1:1
    y = Factor.BTol{i}'*y;
end

y = Factor.SigmaM*y;

for i=1:length(Factor.ATol)
    y = Factor.ATol{i}*y;
end

y = Factor.U*y;

end
