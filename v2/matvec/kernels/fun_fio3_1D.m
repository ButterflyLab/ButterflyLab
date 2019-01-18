function res = fun_fio3_1D(x,k,c)

tmp = 2*pi*( x*k'+(2+0.2*sin(2*pi*x))*(2+0.5*cos(2*pi*k')+c*k'.^2) );
res = complex(cos(tmp),sin(tmp));

end
