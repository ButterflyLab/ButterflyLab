function v = calc_vterm_i(tl,zsrc,jsrc)
% v = calc_vterm_i(tl,zsrc,jsrc)
%
% Part of the tlines calculator (see calc_tlines), given the unit
% shunt current source which belongs to a particular tline, calculate
% the voltage at the left terminal of this tline.
%   tl   - structure with transmission lines parameters and auxiliary
%          data as returned by calc_tlines.
%   zsrc - source coordinate.
%   jsrc - source tline index.

% Possible, but suboptimal way to calculate the voltage
%v = calc_vi(tl.z(jsrc),jsrc,zsrc,jsrc);

% The voltage at terminal is calculated based on the reciprocity
% relations, instead we calculate voltage at the source point due
% to a current source at the terminal. 
% Admittances at the immediate left (Yls) and right (Ygr) points
% of the left terminal of the source tline (terminal jsrc)
Yls=tl.Y0(:,jsrc).*(tl.Gls(:,jsrc)-1)./(tl.Gls(:,jsrc)+1);
Ygr=tl.Y0(:,jsrc-1).*(1-tl.Ggr(:,jsrc-1))./(1+tl.Ggr(:,jsrc-1));
% Find voltage at the terminal due to a unit shunt current source at
% this terminal.
vts=1./(Ygr-Yls);
% vt_vsrc = vsrc/vt, where vsrc is the voltage at the source point
% and vt is the voltage at terminal
ex1 = exp(-2*tl.k(:,jsrc).*(tl.z(:,jsrc+1)-zsrc));
ex2 = exp(-tl.k(:,jsrc).*(zsrc-tl.z(:,jsrc)));
vt_vsrc=(1+tl.Ggr(:,jsrc).*ex1).*ex2./(1+tl.Ggr(:,jsrc).*tl.t(:,jsrc));
% Now, by reciprocity, find the voltage at terminal due to a current
% source at the source point.
v = vts.*vt_vsrc;
