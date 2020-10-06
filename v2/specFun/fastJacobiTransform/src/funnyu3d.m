function M = funnyu3d(rs,cs,n,da,db)
rs = rs*1.000;
cs = cs*1.000;
nt = zeros(nthroot(n,3),1);
M = extrjac3(nt,rs,cs,-1,da,db);
end