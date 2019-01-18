function p = polar_k2p(k)

Dim = size(k,2);

if Dim == 1
    p = k;
    return;
end

% Cartesian coordinates to Polar coordinates
if Dim == 2
    k1 = k(:,1);
    k2 = k(:,2);
    p1 = sqrt(k1.^2+k2.^2);
    p2 = atan2(k2,k1);
    p2(p2<0) = p2(p2<0)+2*pi;
    p = [p1 p2];
    return;
end

% Cartesian coordinates to Spherical coordinates
if Dim == 3
    k1 = k(:,1);
    k2 = k(:,2);
    k3 = k(:,3);
    p1 = sqrt(k1.^2+k2.^2+k3.^2);
    p2 = atan2(sqrt(k1.^2+k2.^2),k3);
    p2(p2<0) = p2(p2<0)+2*pi;
    p3 = atan2(k2,k1);
    p = [p1 p2 p3];
    return;
end

if Dim > 3
    assert( Dim <= 3, 'Polar_k2p does not support data with dimension higher than 3.');
end

end
