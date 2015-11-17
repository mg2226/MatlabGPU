function y=hlpr0_c_vec(l,p,rr,root,aa,R2)
%function y=hlpr0_c_vec(l,p,rr,root,aa,R2)
%
% This is the vectorized version of hlpr0_c.m
% hlpr0_c.m compute radial basis function hlp(r) when inner radius R1 is zero
% modified from /home/rupee/a/doerschu/radial_basis_funcs2/Hlpr0.c
%
% Yili Zheng
% 05/19/2008
%
% Caustion: hlpr0_c_vec.m and hlpr0_c.m need to be used in
% conjunction with jroots_c.m and other *_c.m which are translated
% from their C counterparts.
%
% variable naming convention: R2=Rmax, r=rr, aa=normp
%
% Please note that the _c version  differs from Seunghee's version in two 
%
% 1) Seunghee's roots = _c version's roots/Rmax;
% hlp0_c divides root*rr by R2: 
%   [sj,sjp,sy,syp]=sphbes_vec(l,root(l+1,p).*rr(rr_nz_ind)/R2);
% but Seunghee's hlpr0.m doens't divide root*r by R2 because root has
% already been divided by R2 in jroots.m:
%   [sj,sjp,sy,syp]=sphbes(l,root(l+1,p)*r);
% .
% This difference would not change the final answer.
% 
% 2) aa in Seunghee's code is always positive but not in the _c version.
%
% The C version code works because the signs of the coefficients are the
% opposite in both hlpk0_c.m and hlpr0_c.m for the "l,p" pairs that have 
% negative "df".

rr_nz_ind = (rr>=0 & rr<=R2);
rr_zero_ind = ~rr_nz_ind; %If redundant below, then this is also redundant.
[sj,sjp,sy,syp]=sphbes_vec(l,root(l+1,p).*rr(rr_nz_ind)/R2);
y2=aa(l+1,p).*sj;

Nr=size(rr,1);
y=zeros(Nr,1);

y(rr_zero_ind)=0.0; %Redundant with the initialization by 'zeros'?
y(rr_nz_ind)=y2;
