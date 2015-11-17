function vobj=virusobj_set_2Dreciprocal(vobj,vk)
%function vobj=virusobj_set_2Dreciprocal(vobj,vk)
%vobj: virusobj
%vk(1:Ny,1:2): list of reciprocal-space image frequency vectors
%adds the following fields to vobj:
%vobj{eta}.unique_lp, vobj{eta}.map_lp2unique, vobj{eta}.map_unique2lp, vobj{eta}.Htable

%Set up the table of reciprocal-space radial basis functions.  These
%values are used in the projection-slice theorem so the $\vec k$ of
%interest is a rotation of $(\vec\kappa 0)^T$.  Therefore,
%$|\vec k|=|\vec\kappa|$ and it the the two components of $\vec\kappa$
%that are passed in the argument vk.

if isempty(vk)
  fprintf(1,'virusobj_set_2Dreciprocal: vk==[], no change to vobj\n');
  return;
end

if ndims(vk)~=2
  error('virusobj_set_2Dreciprocal: ndims(vk) %d ~= 2\n',ndims(vk));
end
if size(vk,2)~=2
  error('virusobj_set_2Dreciprocal: size(vk,2) %d ~= 2\n',size(vk,2));
end

fprintf(1,'virusobj_set_2Dreciprocal: vobj{eta}.R1<0.0 means use H_{l,p}(r) on [0,R_2)\n');
fprintf(1,'virusobj_set_2Dreciprocal: vobj{eta}.R1>=0.0 means use H_{l,p}(r) on (R_1,R_2) with R_1>=0\n');

fprintf(1,'virusobj_set_2Dreciprocal: the reciprocal-space radial basis functions are the same for all values of vobj{eta}.BasisFunctionType: \n');
for eta=1:length(vobj)
  fprintf(1,'virusobj_set_2Dreciprocal: vobj{eta=%d}.BasisFunctionType=%d\n',eta,vobj{eta}.BasisFunctionType);
end

for eta=1:length(vobj)

  if vobj{eta}.R1>=0.0
    error('virusobj_set_2Dreciprocal: Only R1<0.0 (not R1>=0.0) is implemented. R1 %g\n',vobj{eta}.R1);
  end
  [vobj{eta}.unique_lp, vobj{eta}.map_lp2unique, vobj{eta}.map_unique2lp] = ...
                                 unique([vobj{eta}.clnp.il vobj{eta}.clnp.ip], 'rows');
  vobj{eta}.Htable=set_H0table_c(vobj{eta}.unique_lp(:,1), vobj{eta}.unique_lp(:,2), ...
                                 vk(:,1), vk(:,2), vobj{eta}.R2);
  lmax=max(vobj{eta}.clnp.il);
  ftable=factorial(0:2*lmax);
  %ftable=zeros(2*lmax+1,1);
  %for ii=1:2*lmax+1
  %    ftable(ii)=factorial(ii-1);
  %end
  vobj{eta}.Ftable=ftable;
  if ~isreal(vobj{eta}.Htable)
    error('virusobj_set_2Dreciprocal: eta %d: Htable is not real!\n',eta);
  end

end
