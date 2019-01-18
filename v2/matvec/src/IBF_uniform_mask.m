function maskFactor = IBF_uniform_mask(fun,xx,kk,fun_mask,r,tol)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.

if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end

Dim = size(xx,2);
xbox = zeros(2,Dim);
kbox = zeros(2,Dim);

N = min(size(xx,1),size(kk,1));

[Nx,~] = size(xx);
[Nk,~] = size(kk);
Nx = Nx^(1/Dim)*ones(1,Dim);
Nk = Nk^(1/Dim)*ones(1,Dim);
xbox(1,:) = min(xx);
xbox(2,:) = max(xx)+(max(xx)-min(xx))./Nx;
kbox(1,:) = min(kk);
kbox(2,:) = max(kk)+(max(kk)-min(kk))./Nk;

% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,r)
levels = max(0,min(floor(log2([Nx Nk])-max(3,ceil(log2(r)))))); 

xxset = cell(levels+1,1);
xidxset = cell(levels+1,1);

xxboxidx = zeros(size(xx));
npxx = 2^levels*ones(1,Dim);
xxset{levels+1} = cell(prod(npxx),1);
xidxset{levels+1} = cell(prod(npxx),1);
for i = 1:Dim
    edges = linspace(xbox(1,i),xbox(2,i),npxx(i)+1);
    [~,xxboxidx(:,i)] = histc(xx(:,i),edges);
end
[xxboxidx,xxidx] = sortrows(xxboxidx,Dim:-1:1);
[xxC,xxIA,~] = unique(xxboxidx,'rows','stable');
xxC = [xxC;zeros(1,size(xxC,2))];
xxIA = [xxIA;size(xx,1)+1];
for itx = 1:prod(npxx)
    x = BF_idx2vec(npxx,itx);
    if any(x ~= xxC(itx,:))
        xxC = [xxC(1:itx-1,:);x;xxC(itx:end,:)];
        xxIA = [xxIA(1:itx-1);xxIA(itx);xxIA(itx:end)];
    end
end
for itx = 1:prod(npxx)
    xxset{levels+1}{itx} = xx(xxidx(xxIA(itx):xxIA(itx+1)-1),:);
    xidxset{levels+1}{itx} = xxidx(xxIA(itx):xxIA(itx+1)-1);
end

for lvl = levels-1:-1:0
    npxx = 2^lvl;
    xxset{lvl+1} = cell(prod(npxx),1);
    xidxset{lvl+1} = cell(prod(npxx),1);
    for itx = 1:prod(npxx)
        xxset{lvl+1}{itx} = [];
        xidxset{lvl+1}{itx} = [];
        for itc = BF_childidx(npxx,itx)
            xxset{lvl+1}{itx} = [ xxset{lvl+1}{itx}; xxset{lvl+2}{itc}];
            xidxset{lvl+1}{itx} = [ xidxset{lvl+1}{itx}; xidxset{lvl+2}{itc}];
        end
    end
end

kkset = cell(levels+1,1);
kidxset = cell(levels+1,1);

kkboxidx = zeros(size(kk));
npkk = 2^levels*ones(1,Dim);
kkset{levels+1} = cell(prod(npkk),1);
kidxset{levels+1} = cell(prod(npkk),1);
for i = 1:Dim
    edges = linspace(kbox(1,i),kbox(2,i),npkk(i)+1);
    [~,kkboxidx(:,i)] = histc(kk(:,i),edges);
end
[kkboxidx,kkidx] = sortrows(kkboxidx,Dim:-1:1);
[kkC,kkIA,~] = unique(kkboxidx,'rows','stable');
kkC = [kkC;zeros(1,size(kkC,2))];
kkIA = [kkIA;size(kk,1)+1];
for itk = 1:prod(npkk)
    k = BF_idx2vec(npkk,itk);
    if any(k ~= kkC(itk,:))
        kkC = [kkC(1:itk-1,:);k;kkC(itk:end,:)];
        kkIA = [kkIA(1:itk-1);kkIA(itk);kkIA(itk:end)];
    end
end
for itk = 1:prod(npkk)
    kkset{levels+1}{itk} = kk(kkidx(kkIA(itk):kkIA(itk+1)-1),:);
    kidxset{levels+1}{itk} = kkidx(kkIA(itk):kkIA(itk+1)-1);
end

for lvl = levels-1:-1:0
    npkk = 2^lvl;
    kkset{lvl+1} = cell(prod(npkk),1);
    kidxset{lvl+1} = cell(prod(npkk),1);
    for itk = 1:prod(npkk)
        kkset{lvl+1}{itk} = [];
        kidxset{lvl+1}{itk} = [];
        for itc = BF_childidx(npkk,itk)
            kkset{lvl+1}{itk} = [ kkset{lvl+1}{itk}; kkset{lvl+2}{itc}];
            kidxset{lvl+1}{itk} = [ kidxset{lvl+1}{itk}; kidxset{lvl+2}{itc}];
        end
    end
end

% Check Oscillatory

queueCell = cell(1);
maskFactor = cell(1);
it = 1;
sp = 1;
ep = 1;
queueCell{1}.lvl = 0;
queueCell{1}.itx = 1;
queueCell{1}.itk = 1;

while sp <= ep
    lvl = queueCell{sp}.lvl;
    itx = queueCell{sp}.itx;
    itk = queueCell{sp}.itk;
    cflag = BF_checkOscillatory(fun_mask, xxset{lvl+1}{itx}, kkset{lvl+1}{itk});
    if cflag == 1
        [Factor,~] = IBF_uniform(fun, ...
            xxset{lvl+1}{itx},kkset{lvl+1}{itk},r,tol,N);
        maskFactor{it}.xidx = xidxset{lvl+1}{itx};
        maskFactor{it}.kidx = kidxset{lvl+1}{itk};
        maskFactor{it}.isfactor = 1;
        maskFactor{it}.factor = Factor;
        it = it+1;
    elseif cflag == 0
        if lvl < levels
            npxx = 2^lvl;
            npkk = 2^lvl;
            for itcx = BF_childidx(npxx,itx)
                for itck = BF_childidx(npkk,itk)
                    ep = ep+1;
                    queueCell{ep}.lvl = lvl+1;
                    queueCell{ep}.itx = itcx;
                    queueCell{ep}.itk = itck;
                end
            end
        else
            maskFactor{it}.xidx = xidxset{lvl+1}{itx};
            maskFactor{it}.kidx = kidxset{lvl+1}{itk};
            maskFactor{it}.isfactor = 0;
            maskFactor{it}.factor = ...
                fun(xxset{lvl+1}{itx},kkset{lvl+1}{itk});
            it = it+1;
        end
    else
        maskFactor{it}.xidx = xidxset{lvl+1}{itx};
        maskFactor{it}.kidx = kidxset{lvl+1}{itk};
        maskFactor{it}.isfactor = 0;
        maskFactor{it}.factor = ...
            fun(xxset{lvl+1}{itx},kkset{lvl+1}{itk});
        it = it+1;
    end
    sp = sp+1;
end

end