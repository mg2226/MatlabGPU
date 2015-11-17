function projimage=test_sphericalobj2(minindex,maxindex,delta,rho0,R2)
%function projimage=test_sphericalobj2(minindex,maxindex,delta,rho0,R2)
%compute the projection of a spherical ball of constant electron scattering intensity rho0 and radius R2
%calculation is done by symbolic formula

if min(size(minindex))~=1 || max(size(minindex))~=2
  error('minindex dimensions');
end
if min(size(maxindex))~=1 || max(size(maxindex))~=2
  error('maxindex dimensions');
end

xrange=[minindex(1):maxindex(1)];
yrange=[minindex(2):maxindex(2)];

[x,y]=ndgrid(xrange,yrange);
rsqu=x.^2+y.^2;
ii=find(rsqu<=(R2/delta)^2);
projimage=zeros(size(x));
projimage(ii)=rho0*2*delta*sqrt((R2/delta)^2 - rsqu(ii));
