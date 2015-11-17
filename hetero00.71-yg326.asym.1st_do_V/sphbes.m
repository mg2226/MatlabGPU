function [sj,sjp,sy,syp]=sphbes(n,x)
%function [sj,sjp,sy,syp]=sphbes(n,x)
%compute spherical bessel function J-Y, and their derivatives
%modified from /home/rupee/a/doerschu/radial_basis_funcs2/sphebes.real.c
%incomplete treatment for x=0 

if n<0 || x<0.0
  error('bad arguments in sphbes : n=%d, x=%g.\n',n,x);
end
if nargout~=4
  error('number of outputs should be two or four.\n');
end

if x==0
  if n==0
    sj=1.0;
  else
    sj=0.0;
  end
  sjp=0.0;
  sy=0.0;
  syp=0.0;

else
  order=n+0.5;
  sj=sqrt(pi/(2*x))*besselj(order,x);
  sy=sqrt(pi/(2*x))*bessely(order,x);
  %by recursion formula
  sjp=(n/x)*sj-sqrt(pi/(2*x))*besselj(order+1,x);
  syp=(n/x)*sy-sqrt(pi/(2*x))*bessely(order+1,x);
end