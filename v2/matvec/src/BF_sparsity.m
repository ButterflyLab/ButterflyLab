function B = BF_sparsity(A)
% This code visualizes the sparsity of the factors stored in A.

B.V = sp(A.V);
B.HTol = cell(1,length(A.HTol));
B.GTol = cell(1,length(A.GTol));
for i=length(A.HTol):-1:1
    B.HTol{i} = sp(A.HTol{i});
end
B.M = sp(A.M);
for i=1:length(A.GTol)
    B.GTol{i} = sp(A.GTol{i});
end
B.U = sp(A.U);

end

function B = sp(A)
B = zeros(size(A));
pos = find(abs(A)>1e-10);
B(pos) = 1;
end