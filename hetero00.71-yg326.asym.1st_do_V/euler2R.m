function R=euler2R(t)
%function R=euler2R(t)
%
% modified from euler2R.c
% M. E. Rose, Elementary Theory of Angular Momentum, p. 50
% and p. 65 (Eq. 4.43). */
%
% Yili Zheng
% 10/31/2007
%
%  /* t=(alpha,beta,gamma) */

ca = cos(t(1));
sa = sin(t(1));
cb = cos(t(2));
sb = sin(t(2));
cg = cos(t(3));
sg = sin(t(3));

R = zeros(3,3);

R(1,1) = ca * cb * cg - sa * sg;
R(1,2) = sa * cb * cg + ca * sg;
R(1,3) = -sb * cg;

R(2,1) = -ca * cb * sg - sa * cg;
R(2,2) = -sa * cb * sg + ca * cg;
R(2,3) = sb * sg;

R(3,1) = ca * sb;
R(3,2) = sa * sb;
R(3,3) = cb;
