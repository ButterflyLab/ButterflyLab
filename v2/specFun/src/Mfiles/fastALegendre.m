function [valp,valq,alpha,alphader,vallogp,vallogq] = fastALegendre(degree,order,t)
% !  This file contains code for evaluating the scaled, normalized associated Legendre
% !  functions of the first and seconds kinds
% !
% !
% !     (  (degree + 1/2) Gamma(degree+order+1) )
% !     (  ------------------------------------ )^(1/2)  P_degree^{-order}(\cos(t)) \sqrt{\sin(t)}       (1)
% !     (        Gamma(degree-order+1)          )
% !
% !  and
% !
% !     2  (  (degree + 1/2) Gamma(degree+order+1) )
% !   ---- (  ------------------------------------ )^(1/2)  Q_degree^{-order}(\cos(t)) \sqrt{\sin(t)}    (2)
% !    Pi  (        Gamma(degree-order+1)          )
% !
% !  on the interval (0,\pi).  It runs in time independent of the degree degree and order order
% !  and can be applied when
% !
% !     0 <= degree <= 1,000,000 and 0 <= order <= degree.                                          (3)
% !
% !  There is also code for calculating the roots of (1) or (2) on the interval
% !  (0,\pi) which can be applied when
% !
% !     10 <= degree <= 1,000,000 and 0 <= order <= degree.
% !
% !  The time taken to evaluate each root is independent of degree and order.
% !
% !  For the most part, these tasks are accomplished using precomputed tables stored in
% !  the files alegendre_data.bin1 and alegendre_data.bin2.   These tables must be read into
% !  memory by the routine  "alegendre_eval_init" before any other subroutine is called.
% !  The algorithm used for the evaluation of (1) and (2) is described in the preprint
% !
% !    James Bremer, "An algorithm for the numerical evaluation of the associated Legendre
% !    functions in time independent of degree and order."  arXiv:1707.03287
% !
% !  and the algorithm used for the calculation of roots will be described in an upcoming
% !  article.
% !
% !  The principal subroutine for the evaluation of (1) and (2) is called  alegendre_eval.
% !  We call the subset of R^3
% !
% !    { (degree,order,t) : degree >= 0, 0 <= order <= 1/2 and 0 <t <= pi }  \cup
% !    { (degree,order,t) : degree >= 0, order > 1/2 and  tp <= t <= \pi/2 }  \cup                    (4)
% !    { (degree,order,t) : degree >0 order > 1/2 and pi/2 < t <= pi - tp  },
% !
% !  where tp the turning point tp  = arcsin( sqrt(order^2-1/4) / (degree+1/2) ), the
% !  oscillatory region.  When (degree,order,t) is in the oscillatory region, in addition
% !  to the values of (1) and (2), alegendre_eval returns the values of a nonoscillatory
% !  phase function alpha for  the associated Legendre differential equation and its
% !  derivative.
% !
% !  The nonoscillatory region is
% !
% !    { (degree,order,t) : degree >=0, order > 1/2 and 0 < t < tp }  \cup                            (5)
% !    { (degree,order,t) : degree >=0, order > 1/2 and pi - tp < t  <pi}
% !

% check inputs
%if (t<=0 | t>=pi), error('t out of range'); end;
%if (degree<0 | degree > 1000000), error('degree out of range'); end;
%if (order<0 | order > degree), error('order out of range'); end;
alpha = zeros(numel(t),numel(degree));
alphader = zeros(numel(t),numel(degree));
vallogp = zeros(numel(t),numel(degree));
vallogq = zeros(numel(t),numel(degree));
valp = zeros(numel(t),numel(degree));
valq = zeros(numel(t),numel(degree));
t = t(:);
for cnt = 1:numel(degree)
    [alpha(:,cnt),alphader(:,cnt),vallogp(:,cnt),vallogq(:,cnt),valp(:,cnt),valq(:,cnt)] = fastALegendre_mex(degree(cnt),order,t);
end
end
