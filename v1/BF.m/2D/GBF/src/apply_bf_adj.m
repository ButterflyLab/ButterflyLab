function y = apply_bf_adj(Factor, y)

if( size(y,2)==0 )
    return;
end
y = Factor.U'*y;

for i=length(Factor.ATol):-1:1
    y = Factor.ATol{i}'*y;
end

y = Factor.SigmaM'*y;

for i=1:length(Factor.BTol)
    y = Factor.BTol{i}*y;
end

y = Factor.V*y;

end
