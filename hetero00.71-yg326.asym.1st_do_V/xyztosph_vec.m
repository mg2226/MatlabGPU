function vecsph= xyztosph_vec(x)
%vecsph=[r, theta, phi];
%function [r,theta,phi]= xyztosph_vec(x);
%convert cartesian coordinates to spherical coordinates
%vectorized version of xyztosph.m, which is modified from xyztosph.c

if size(x,2)~=3
   error('x should be a n-by-3 array.');
end

r=zeros(size(x,1),1);
theta=zeros(size(x,1),1);
phi=zeros(size(x,1),1);

r(:) = sqrt(x(:,1).*x(:,1)+x(:,2).*x(:,2)+x(:,3).*x(:,3));

for i=1:size(x,1)
  if r(i)==0.0
    theta(i)=0.0;
    phi(i)=0.0;
  else
    theta(i) = acos(x(i,3)/r(i));
    if x(i,1)==0.0 && x(i,2)==0.0
      phi(i) = 0.0;
    else
      phi(i) = atan2(x(i,2), x(i,1));
    end
  end
end

vecsph=[r,theta,phi];
