function vobj=virusobj_read(fn_clnp,fn_nu,fn_q)
%function vobj=virusobj_read(fn_clnp,fn_nu,fn_q)
%fn_clnp{1:Neta} containing filenames as strings -- mean of weights
%fn_nu{1:Neta} containing filenames as strings -- diagonal (as a vector) of covariance of weights
%fn_q{1:Neta} containing filenames as strings -- prior probabilities, one per file

Neta=length(fn_clnp);

if length(fn_nu)~=Neta | length(fn_q)~=Neta
  error('virusobj: one or more of fn_nu and fn_q are not Neta %d long\n',Neta);
end

if Neta==0
  fprintf(1,'virusobj_read: Neta = 0\n');
  vobj=[];
  return;
end

if ~iscellstr(fn_clnp) | ~iscellstr(fn_nu) | ~iscellstr(fn_q)
  error('virusobj: one or more of fn_clnp, fn_nu, and fn_q are not a cell array of strings\n');
end

%allocate vobj:
vobj=cell(Neta, 1);

for eta=1:Neta

  %load file names:
  vobj{eta}.clnp_fn=fn_clnp{eta};
  vobj{eta}.nu_fn=fn_nu{eta}; %PCD 2013-06-18 addition to YZ structure
  vobj{eta}.q_fn=fn_q{eta}; %PCD 2013-06-18 addition to YZ structure

  %read vobj{eta}.clnp_fn:
  [header1,header2,vobj{eta}.clnp]=rd_clnp(vobj{eta}.clnp_fn);
  %the previous command set the following:
  %vobj{eta}.clnp.il
  %vobj{eta}.clnp.in
  %vobj{eta}.clnp.ip
  %vobj{eta}.clnp.optflag
  %vobj{eta}.clnp.c
  vobj{eta}.cbar=vobj{eta}.clnp.c; %for convenience
  vobj{eta}.BasisFunctionType=header1;
  vobj{eta}.R1=header2(1);  % header2(1)=R1, header2(2)=R2;
  vobj{eta}.R2=header2(2);  % header2(1)=R1, header2(2)=R2;
  fprintf(1,'virusobj_read: eta %d, R1 %g, R2 %g\n',eta,header2(1),header2(2));
  %test to make sure that there are no repeated weights
  indices=[vobj{eta}.clnp.il vobj{eta}.clnp.in vobj{eta}.clnp.ip];
  uniqueindices=unique(indices,'rows');
  if size(indices,1)~=size(uniqueindices,1)
    error('virusobj_read: repeated weights\n');
  end
  clear indices uniqueindices

  %read vobj{eta}.nu_fn:
  [vobj{eta}.nu,count]=rd_nu(vobj{eta}.nu_fn);
  if count~=0 && count~=length(vobj{eta}.cbar)
    error('virusobj_read: lengths of nu and cbar are different\n');
  end

  %read vobj{eta}.q_fn:
  [vobj{eta}.q,count]=rd_q(vobj{eta}.q_fn);

end

tmp=0;
for eta=1:Neta
  tmp=tmp+vobj{eta}.q;
end
if abs(tmp-1.0)>1.0e-12
  error('virusobj_read: vobj{eta}.q summed over eta does not equal 1\n');
end
