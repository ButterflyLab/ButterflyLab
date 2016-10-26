function k = polar_p2k(p)

Dim = size(p,2);

if Dim == 1
    k = p;
    return;
end

% Polar coordinates to Cartesian coordinates
if Dim == 2
    p1 = p(:,1);
    p2 = p(:,2);
    
    k1 = p1.*cos(p2);
    k2 = p1.*sin(p2);
    k = [k1 k2];
    return;
end

%  Spherical coordinates to Cartesian coordinates
if Dim == 3
    p1 = p(:,1);
    p2 = p(:,2);
    p3 = p(:,3);
    
    k1 = p1.*sin(p2).*cos(p3);
    k2 = p1.*sin(p2).*sin(p3);
    k3 = p1.*cos(p2);
    k = [k1 k2 k3];
    return;
end

if Dim > 3
    assert( Dim <= 3, 'Polar_p2k does not support data with dimension higher than 3.');
end

end
