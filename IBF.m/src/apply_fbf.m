function y = apply_fbf(Factor, y)

y = Factor.V*y;

for i=length(Factor.HTol):-1:1
    y = Factor.HTol{i}*y;
end

y = Factor.M*y;

for i=1:length(Factor.GTol)
    y = Factor.GTol{i}*y;
end

y = Factor.U*y;

end
