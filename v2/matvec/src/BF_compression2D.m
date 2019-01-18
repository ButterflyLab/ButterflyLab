function [WPreSpr,CPreSpr] = BF_compression2D(WPreSpr, CPreSpr, W, kk, kkgid, kbox, r, tol, level, NL, XYmesh)

H = size(W{1,1},1);
[km1,km2] = size(W);


if( H <= 0 || level > NL )
    if( level <= NL && H <= 2*r )
        COffset = CPreSpr(level).Offset;
        CHeight = CPreSpr(level).Height;
        CWidth = CPreSpr(level).Width;
        totalW = 0;
        for k1=1:km1
            for k2=1:km2
                totalW = totalW+size(W{k1,k2},2);
            end
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
        for k1=1:km1
            for k2=1:km2
                [X,Y] = meshgrid(kkgid,WWidth+(1:size(W{k1,k2},2)));
                X = X';
                Y = Y';
                idx = WOffset+(1:numel(X));
                WPreSpr.XT(idx) = X(:);
                WPreSpr.YT(idx) = Y(:);
                WPreSpr.ST(idx) = W{k1,k2}(:);
                WWidth = WWidth+size(W{k1,k2},2);
                if(~isempty(idx))
                    WOffset = idx(end);
                end
                
            end
        end
        WPreSpr.Offset = WOffset;
        WPreSpr.Width = WWidth;
        return;
    end
    [WPreSpr,CPreSpr] = compression2D(WPreSpr, CPreSpr, W, kk, kkgid, kbox, r, tol, level+1, NL, XYmesh);
    return;
end

WW = cell(2,2);
WW{1,1} = cell(km1/2,km2/2);
WW{1,2} = cell(km1/2,km2/2);
WW{2,1} = cell(km1/2,km2/2);
WW{2,2} = cell(km1/2,km2/2);
CS = cell(2,2,km1,km2);

kidx = bf_prep(kk,kbox,2,2);
for k1=1:2:km1
    for k2=1:2:km2
        hh = size(W{k1,k2},1);
        r1 = size(W{k1,k2},2);
        r2 = size(W{k1,k2+1},2);
        r3 = size(W{k1+1,k2},2);
        r4 = size(W{k1+1,k2+1},2);
        for i=1:2
            for j=1:2
                Id = kidx{i,j};
                if(hh>0)
                    chunk = [W{k1,k2}(Id,:), W{k1,k2+1}(Id,:), ...
                        W{k1+1,k2}(Id,:), W{k1+1,k2+1}(Id,:) ];
                    [Wtmp,S,V] = svdrand(chunk,r,tol);
                    WW{i,j}{(k1+1)/2,(k2+1)/2} = Wtmp*S;
                    CS{i,j,k1,k2} = V(1:r1,:)';
                    CS{i,j,k1,k2+1} = V(r1+(1:r2),:)';
                    CS{i,j,k1+1,k2} = V(r1+r2+(1:r3),:)';
                    CS{i,j,k1+1,k2+1} = V(r1+r2+r3+(1:r4),:)';
                end
            end
        end
    end
end

COffset = CPreSpr(level).Offset;
CHeight = CPreSpr(level).Height;
CWidth = CPreSpr(level).Width;
totalW = CWidth;
currentW = zeros(km1,km2);
for k1=1:km1
    for k2=1:km2
        currentW(k1,k2) = totalW;
        totalW = totalW+size(W{k1,k2},2);
    end
end

totalH = CHeight;
offset = COffset;
for z1=1:2
    for z2=1:2
        for kk1=1:2:km1
            for kk2=1:2:km2
                for k1 = kk1:kk1+1
                    for k2 = kk2:kk2+1
                        tmpM = CS{z1,z2,k1,k2};
                        if(min(size(tmpM))>0)
                            X = XYmesh{size(tmpM,1),size(tmpM,2),1}+totalH;
                            Y = XYmesh{size(tmpM,1),size(tmpM,2),2}+currentW(k1,k2);
                            idx = offset+1:offset+numel(X);
                            CPreSpr(level).XT(idx) = X(:);
                            CPreSpr(level).YT(idx) = Y(:);
                            CPreSpr(level).ST(idx) = tmpM(:);
                            if(~isempty(idx))
                                offset = idx(end);
                            end
                        end
                    end
                end
                totalH = totalH + size(CS{z1,z2,kk1,kk2},1);
            end
        end
    end
end

CPreSpr(level).Offset = offset;
CPreSpr(level).Height = totalH;
CPreSpr(level).Width = totalW;

k1os = kbox(1,1);
k1len = kbox(1,2)-kbox(1,1);
k2os = kbox(2,1);
k2len = kbox(2,2)-kbox(2,1);
for z1=1:2
    for z2=1:2
        Id = kidx{z1,z2};
        ksubbox = [ k1os+(z1-1)*k1len/2, k1os+z1*k1len/2; ...
            k2os+(z2-1)*k2len/2, k2os+z2*k2len/2];
        [WPreSpr,CPreSpr] = compression2D(WPreSpr, CPreSpr, WW{z1,z2},...
            kk(Id,:), kkgid(Id), ksubbox, r, tol, level+1, NL, XYmesh);
    end
end

end
