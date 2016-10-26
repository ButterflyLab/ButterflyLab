function res = funFT(x,k)

tmp = (2*pi)* (x*k');
res = complex(cos(tmp),-sin(tmp));

end
