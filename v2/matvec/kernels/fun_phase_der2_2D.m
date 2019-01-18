function tmp = fun_phase_der2_2D(x,k)

xk = x(:,2)*ones(1,size(k,1));
sx = (2 + sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)))/32;
cx = (2 + cos(2*pi*x(:,1)).*cos(2*pi*x(:,2)))/32;
kr = (1./sqrt(sx.^2*(k(:,1).^2)' + cx.^2*(k(:,2).^2)')/2).*( 2*cx.^2*(k(:,2))' );

tmp = (2*pi)* (xk + kr);

end
