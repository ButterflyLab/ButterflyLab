function ks = pbf_p2k(N,ps)

% Polar coordinates to Cartesian coordinates
p1 = ps(:,1);
p2 = ps(:,2);

k1 = sqrt(2)/2*N* p1.*cos(2*pi*p2);
k2 = sqrt(2)/2*N* p1.*sin(2*pi*p2);
ks = [k1 k2];

end
