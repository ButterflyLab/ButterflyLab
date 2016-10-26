function [WPreSpr, CPreSpr] = compression1D(WPreSpr, CPreSpr, W, kk, kkgid, kbox, r, tol, level, NL)

H = size(W{1},1);
[km,~] = size(W);

if( H <= 2*r || level > NL )
    if( level <= NL && H <= 2*r )
        COffset = CPreSpr(level).Offset;
        CHeight = CPreSpr(level).Height;
        CWidth = CPreSpr(level).Width;
        totalW = 0;
        for k=1:km
            totalW = totalW+size(W{k},2);
        end
        CPreSpr(level).XT(COffset+(1:totalW)) = CHeight+(1:totalW);
        CPreSpr(level).YT(COffset+(1:totalW)) = CWidth+(1:totalW);
        CPreSpr(level).ST(COffset+(1:totalW)) = ones(totalW,1);
        CPreSpr(level).Offset = COffset+totalW;
        CPreSpr(level).Height = CHeight+totalW;
        CPreSpr(level).Width = CWidth+totalW;
    end
    if( level > NL )
        WOffset = WPreSpr.Offset;
        WWidth = WPreSpr.Width;
        for k=1:km
            [X,Y] = meshgrid(kkgid,WWidth+(1:size(W{k},2)));
            X = X';
            Y = Y';
            idx = WOffset+(1:numel(X));
            WPreSpr.XT(idx) = X(:);
            WPreSpr.YT(idx) = Y(:);
            WPreSpr.ST(idx) = W{k}(:);
            WWidth = WWidth+size(W{k},2);
            if(~isempty(idx))
                WOffset = idx(end);
            end
        end
        WPreSpr.Offset = WOffset;
        WPreSpr.Width = WWidth;
        return;
    end
    [WPreSpr,CPreSpr] = compression1D(WPreSpr, CPreSpr, W, kk, kkgid, kbox, r, tol, level+1, NL);
    return;
end

WW = cell(2,1);
WW{1} = cell(km/2,1);
WW{2} = cell(km/2,1);
CS = cell(2,km);

kidx = bf_prep(kk,kbox,2);
for k=1:2:km
    r1 = size(W{k},2);
    r2 = size(W{k+1},2);
    for i=1:2
        Id = kidx{i};
        chunk = [W{k}(Id,:), W{k+1}(Id,:) ];
        [Wtmp,S,V] = svdtrunc(chunk,r,tol);
        WW{i}{(k+1)/2} = Wtmp*S;
        CS{i,k} = V(1:r1,:)';
        CS{i,k+1} = V(r1+(1:r2),:)';
    end
end

COffset = CPreSpr(level).Offset;
CHeight = CPreSpr(level).Height;
CWidth = CPreSpr(level).Width;
totalW = CWidth;
currentW = zeros(km,1);
for k=1:km
    currentW(k) = totalW;
    totalW = totalW+size(W{k},2);
end

totalH = CHeight;
offset = COffset;
for z=1:2
    for kk1=1:2:km
        for k = kk1:kk1+1
            tmpM = CS{z,k};
            [X,Y] = meshgrid(totalH+(1:size(tmpM,1)), ...
                currentW(k)+(1:size(tmpM,2)));
            X = X';
            Y = Y';
            idx = offset+1:offset+numel(X);
            CPreSpr(level).XT(idx) = X(:);
            CPreSpr(level).YT(idx) = Y(:);
            CPreSpr(level).ST(idx) = tmpM(:);
            if(~isempty(idx))
                offset = idx(end);
            end
        end
        totalH = totalH + size(CS{z,kk1},1);
    end
end

CPreSpr(level).Offset = offset;
CPreSpr(level).Height = totalH;
CPreSpr(level).Width = totalW;

kos = kbox(1);
klen = kbox(2)-kbox(1);
for z=1:2
    Id = kidx{z};
    ksubbox = [ kos+(z-1)*klen/2, kos+z*klen/2 ];
    [WPreSpr,CPreSpr] = compression1D(WPreSpr, CPreSpr, WW{z}, kk(Id), kkgid(Id), ksubbox, r, tol, level+1, NL);
end

end
