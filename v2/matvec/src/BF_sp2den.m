function y = BF_sp2den(Factor)
% This code transfers a BF to a dense matrix.

y = full(Factor.V);

for i=length(Factor.HTol):-1:1
    y = full(Factor.HTol{i})*y;
end

y = full(Factor.M)*y;

for i=1:length(Factor.GTol)
    y = full(Factor.GTol{i})*y;
end

y = full(Factor.U)*y;

end
