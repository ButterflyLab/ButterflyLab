function [grid,Lgrid] = fbf_grid(x,npx,chey_grid,box)

Dim = size(x,2);
len = (box(2,:)-box(1,:))./npx;
offset = box(1,:)+(x-1).*len;
Lgrid = chey_grid*len+ones(size(chey_grid))*offset;

I=cell(Dim,1);
for di=1:Dim
    I{di}=Lgrid(:,di);
end
[I{1:Dim}]=ndgrid(I{:});

grid = zeros(length(chey_grid)^Dim,Dim);
for di=1:Dim
    grid(:,di)=I{di}(:);
end

end