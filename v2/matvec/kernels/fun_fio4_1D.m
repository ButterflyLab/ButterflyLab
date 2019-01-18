function res = fun_fio4_1D(x,k,c)
% nonlinearity increases as the range of k increases

tmp = 2*pi*( (x+ c*sin(2*pi*x))*(k'+c*cos(2*pi*k')) );
res = complex(cos(tmp),sin(tmp));

end
