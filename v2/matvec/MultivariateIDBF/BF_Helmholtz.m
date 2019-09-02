% Generation of the grid points for the oscillatory part of the Green's 
% function of a Helmholtz equation.

function [cart, sph, elem] = BF_Helmholtz(r, refinement)
surface = spheresurface(r);
[node, elem] = surface.initmesh();
for i = 1:1:refinement
    [node, elem] = BF_uniformrefine(node, elem);
    node = surface.project(node);
end
[azimuth,elevation,~] = cart2sph(node(:,1),node(:,2),node(:,3));
cart = node;
sph = [azimuth,elevation];
end

function surfacedata = spheresurface(r)
surfacedata = struct('phi', @(p)phi(p,r), 'gradient', @(p)gradient(p,r),...
    'unitoutnormal', @(p)gradient(p,r), 'project', @(p)project(p,r),...
    'Hessian', @(p)Hessian(p,r), 'tanproperator', @(p)tanproperator(p,r),...
    'initmesh',@()initmesh(r));
end

function [node,elem] = initmesh(r)
%Produce a icosahedron mesh.
t=(sqrt(5)-1)/2;
node = [0, 1, t; 0, 1,-t; 1, t, 0; 1,-t, 0; 0,-1,-t; 0,-1, t; ...
        t, 0, 1;-t, 0, 1; t, 0,-1;-t, 0,-1;-1, t, 0;-1,-t, 0];
elem = [6, 2, 0; 3, 2, 6; 5, 3, 6; 5, 6, 7; 6, 0, 7; ...
        3, 8, 2; 2, 8, 1; 2, 1, 0; 0, 1,10; 1, 9,10; ...
        8, 9, 1; 4, 8, 3; 4, 3, 5; 4, 5,11; 7,10,11; ...
        0,10, 7; 4,11, 9; 8, 4, 9; 5, 7,11;10, 9,11] + 1;
node = project(node, r);
end

function p = project(p, r)
% projection function 
L = sqrt(sum(p.^2, 2));
p = p./(L*ones(1,3));
p = r*p;
end

function z = phi(p, r)
% level set function \phi
z = p(:,1).^2+p(:,2).^2+p(:,3).^2-r^2;
end

function n = gradient(p, r)
% the gradient function of \phi
L = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
n = p ./(L*ones(1,3));
end

function H = Hessian(p, r)
% Hessian matrix
H = zeros(3,3,size(p,1));
L = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
L3= L.^3;
d11 = 1./L-p(:,1).^2./L3;
d12 = -p(:,1).*p(:,2)./L3;
d13 = -p(:,1).*p(:,3)./L3;
d22 = 1./L-p(:,2).^2./L3;
d23 = -p(:,2).*p(:,3)./L3;
d33 = 1./L-p(:,3).^2./L3;
for i = 1:size(p,1)
    H(:,:,i)=[d11(i),d12(i),d13(i);d12(i),d22(i),d23(i);...
        d13(i),d23(i),d33(i)];
end
end

function z = tanproperator(p, r)
z = zeros(3,3,size(p,1));
n = gradient(p, r);
for i = 1:size(p,1)
    z(:,:,i) = n(i,:)'*n(i,:); 
end
end