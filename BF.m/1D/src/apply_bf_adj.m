function y = apply_bf_adj(Factor, p)

y = Factor.U'*p;

for i=length(Factor.ATol):-1:1
    y = Factor.ATol{i}'*y;
end

y = Factor.SigmaM'*y;

for i=1:length(Factor.BTol)
    y = Factor.BTol{i}*y;
end

y = Factor.V*y;

end
