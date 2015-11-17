function ht=set_H0table_c(ll,pp,kappa1,kappa2,rmax)
%function ht=set_H0table_c(ll,pp,kappa1,kappa2,rmax)
%computes the table of Hlpk0 : (mkappa', mc') when inner radius is zero.
%modified from set_Htable.m

kappa=sqrt(kappa1.^2+kappa2.^2); %unit of kappa = 1 / (unit of r)
lmax=max(ll);
pmax=max(pp);

if length(ll) ~= length(pp)
  error('The length of ll %d does not equal to the length of pp %d!\n',...
    length(ll), length(pp));
end

%initialize variables for hlpr/hlpk
rtmax=10000.0;
if (rmax > rtmax)
  error('rmax %g > rtmax %g in init_hlp0_Hlp0!\n', rmax, rtmax)
end
[root,normp]=init_hlp0_Hlp0_c(rmax,lmax,pmax,rtmax);

ht=zeros(length(kappa1),length(ll));
for jj=1:length(ll)
  l=ll(jj);
  p=pp(jj);%fprintf(1,'l=%d,p=%d\n',l,p);
  ht(:,jj)=hlpk0_c_vec(l,p,kappa(:),root,rmax);
end