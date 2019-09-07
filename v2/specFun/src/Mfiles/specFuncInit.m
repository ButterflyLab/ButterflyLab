function specFuncInit()
degree = 10; order = 0; rt = 0.5; j=1;
[alphas,alphaders,vallogjs,vallogys,valjs,valys] = fastBessel_mex(degree,rt);
[alpha,alphader,vallogp,vallogq,valp,valq] = fastALegendre_mex(degree,order,rt);
t=alegendreQRoot_mex(degree,order,j);
t=alegendrePRoot_mex(degree,order,j);
[a,alpha,alphader,alphader2]=alegendreTp_mex(degree,order);
[n,m]=alegendreNRoots_mex(degree,order);
t=alegendreInverse_mex(degree,order,2*pi);
[t,w]=alegendreJacobi_mex(order,degree,j);
[alphader]=alegendrePhaseDer_mex(degree,order,rt);
end