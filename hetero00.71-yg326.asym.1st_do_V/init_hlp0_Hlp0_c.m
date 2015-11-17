function [root,normp]=init_hlp0_Hlp0_c(rmax,lmax,pmax,rtmax)
%function [root,normp]=init_hlp0_Hlp0_c(lmax,pmax,rtmax) 
%initialize variables for hlpr0/Hlpk0 

%initialize 
root=zeros(lmax,pmax);
normp=zeros(lmax,pmax);
tmp = sqrt(rmax) * rmax;

for l=0:lmax
  [root(l+1,:), normp(l+1,:)]=jroots_c(l,pmax,rtmax);
  for p = 1:pmax
    normp(l+1,p) = 1.0 / (normp(l+1,p) * tmp); % normp=aa
  end
end

