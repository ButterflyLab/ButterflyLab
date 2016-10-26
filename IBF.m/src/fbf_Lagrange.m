function res = fbf_Lagrange(gs,ts)

[NG,Dim] = size(gs);
NT = size(ts,1);
res = ones(NG^Dim,NT);
for ia=1:NG^Dim
    a = idx2vec(NG*ones(1,Dim),ia);
    gud = cell(Dim,1);
    for di = 1:Dim
        gud{di} = [1:a(di)-1 a(di)+1:NG];
    end
    for b=1:NT
        for di = 1:Dim
            cur = (ts(b,di)-gs(gud{di},di))./(gs(a(di),di)-gs(gud{di},di));
            res(ia,b) = res(ia,b)*prod(cur);
        end
    end
end

end
