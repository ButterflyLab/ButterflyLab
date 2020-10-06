function M = funnyu1d(rs,cs,n,da,db,ts,nu)
rs = rs*1.000;
cs = cs*1.000;
nt = zeros(n,1);
M = extrjac1(nt,rs,cs,-1,da,db,ts,nu);
end
