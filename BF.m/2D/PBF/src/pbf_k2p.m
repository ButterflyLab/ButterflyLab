function ps = pbf_k2p(N,ks)

% Cartesian coordinates to Polar coordinates
k1 = ks(:,1);
k2 = ks(:,2);
kr = sqrt(k1.^2+k2.^2);
p1 = kr / (sqrt(2)/2*N);
ka = atan2(k2,k1);
ka(ka<0) = ka(ka<0)+2*pi;
p2 = ka/(2*pi);
ps = [p1 p2];

end
