function flag = BF_checkOscillatory(fun_mask, xx, kk)

tflag = false;
fflag = false;

Dim = size(xx,2);
flagxx = zeros([2 size(xx)]);
flagkk = zeros([2 size(kk)]);
for it = 1:Dim
    flagxx(1,:,it) = xx(:,it) == min(xx(:,it));
    flagxx(2,:,it) = xx(:,it) == max(xx(:,it));
    flagkk(1,:,it) = kk(:,it) == min(kk(:,it));
    flagkk(2,:,it) = kk(:,it) == max(kk(:,it));
end

xdi = zeros(2^Dim,1);
for di = 1:2^Dim
    biv = BF_de2bi(di-1,Dim)+1;
    vec = ones(size(xx,1),1);
    for it = 1:Dim
        vec = vec & reshape(flagxx(biv(it),:,it),size(xx,1),1);
    end
    xdi(di) = find(vec,1);
end

kdi = zeros(2^Dim,1);
for di = 1:2^Dim
    biv = BF_de2bi(di-1,Dim)+1;
    vec = ones(size(kk,1),1);
    for it = 1:Dim
        vec = vec & reshape(flagkk(biv(it),:,it),size(kk,1),1);
    end
    kdi(di) = find(vec,1);
end

for itx = 1:2^Dim
    for itk = 1:2^Dim
        if ~fun_mask(xx(xdi(itx),:),kk(kdi(itk),:))
            fflag = true;
        else
            tflag = true;
        end
    end
end

if tflag & fflag
    flag = 0;
elseif fflag
    flag = -1;
else
    flag = 1;
end

end