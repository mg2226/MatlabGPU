function projimage=test_sphericalobj1(minindex,maxindex,delta,rho0,R2)
%function projimage=test_sphericalobj1(minindex,maxindex,delta,rho0,R2)
%compute the projection of a spherical ball of constant electron scattering intensity rho0 and radius R2
%calculation is done by creating the object in R^3 and summing over one of the dimensions.

if min(size(minindex))~=1 || max(size(minindex))~=3
  error('minindex dimensions');
end
if min(size(maxindex))~=1 || max(size(maxindex))~=3
  error('maxindex dimensions');
end

xrange=[minindex(1):maxindex(1)];
yrange=[minindex(2):maxindex(2)];
zrange=[minindex(3):maxindex(3)];

[x,y,z]=ndgrid(xrange,yrange,zrange);
rsqu=x.^2+y.^2+z.^2;
ii=find(rsqu<=(R2/delta)^2);
rho=zeros(size(x));
rho(ii)=rho0;
projimage=delta.*sum(rho,3); %does not matter which index is summed but summing on third dimension gives a 2-D result with no index whose range is just 1
