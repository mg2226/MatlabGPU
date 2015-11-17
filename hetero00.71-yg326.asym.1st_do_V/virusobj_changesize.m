function new_vobj=virusobj_changesize(vlmax,vpmax,old_vobj)
%function new_vobj=virusobj_changesize(vlmax,vpmax,old_vobj)
%vlmax(1:Neta): new values for lmax, non-negative integers
%vpmax(1:Neta): new values for pmax, positive integers
%old_vobj{1:Neta}: old virus structures

Neta=length(vlmax);

for eta=1:Neta
  if old_vobj{eta}.BasisFunctionType~=0 && old_vobj{eta}.BasisFunctionType~=1
    error('virusobj_changesize: BasisFunctionType eta %d\n',eta);
  end
end

if Neta~=length(vpmax) | Neta~=length(old_vobj)
  error('virusobj_changesize: not all arguments have the same number of elements\n');
end
if ~isnumeric(vlmax) | ~isnumeric(vpmax)
  error('virusobj_changesize: at least one of vlmax and vpmax is not a numeric array\n');
end
if any(vlmax<0)
  error('virusobj_changesize: at least one negative value in vlmax\n');
end
if any(vpmax<=0)
  error('virusobj_changesize: at least one non-positive value in vlmax\n');
end
if ~iscell(old_vobj)
  error('virusobj: old_vobj is not a cell array\n');
end

%For harmonic angular basis functions, it is not possible to have more than 2l+1 basis functions for each value of l.  Therefore, the maximum number of angular basis functions to order l_{\rm max} is \sum_{l=0}^{l_{\rm max}} (2l+1) = (l_{\rm max} +1)^2.

new_vobj=old_vobj;

for eta=1:Neta

  lmax=vlmax(eta);
  pmax=vpmax(eta);

  llist=zeros((lmax+1)^2,1);
  nlist=zeros((lmax+1)^2,1);

  %the following code causes n to change more rapidly than l
  ln=0;
  for l=0:lmax
    for n=0:2*l
      %selection rule for icosahedrally-symmetric harmonic angular basis functions
      if (old_vobj{eta}.BasisFunctionType==0) || (old_vobj{eta}.BasisFunctionType==1 && slct(l,n)==1)
        ln=ln+1;
        llist(ln)=l;
        nlist(ln)=n;
      end
    end
  end

  fprintf(1,'virusobj_changesize: eta %d ln %d\n',eta,ln);

  new_il=zeros(ln*pmax,1);
  new_in=zeros(ln*pmax,1);
  new_ip=zeros(ln*pmax,1);
  new_optflag=ones(ln*pmax,1); %default for optflag is 1
  new_c=zeros(ln*pmax,1); %default for clnp is 0.0
  if ~isempty(old_vobj{eta}.nu)
    new_nu=zeros(ln*pmax,1); %default for nu is 0.0
  else
    new_nu=[];
  end

  %the following code causes p to change most rapidly
  lnp=0;
  for ii=1:ln
    for p=1:pmax
      lnp=lnp+1;
      new_il(lnp)=llist(ii);
      new_in(lnp)=nlist(ii);
      new_ip(lnp)=p;
    end
  end

  old_il=old_vobj{eta}.clnp.il;
  old_in=old_vobj{eta}.clnp.in;
  old_ip=old_vobj{eta}.clnp.ip;
  old_optflag=old_vobj{eta}.clnp.optflag;
  old_c=old_vobj{eta}.clnp.c;
  old_nu=old_vobj{eta}.nu;

  for ii=1:lnp
%    fprintf(1,'virusobj_changesize: ii %d seeking (%d %d %d)\n',ii,new_il(ii),new_in(ii),new_ip(ii));
    logical_l=(old_il==new_il(ii));
    logical_n=(old_in==new_in(ii));
    logical_p=(old_ip==new_ip(ii));
    logical_lnp=(logical_l & logical_n & logical_p);
    jj=find(logical_lnp==1);
    if length(jj)>1 
      error('virusobj_changesize: ii %d found duplicates\n',ii);
    end
    if length(jj)==1
      new_optflag(ii)=old_optflag(jj);
      new_c(ii)=old_c(jj);
      if ~isempty(old_vobj{eta}.nu)
        new_nu(ii)=old_nu(jj);
      end
    end
    %ok if no corresponding old value was found, default is already set
  end

  %replace existing values
  new_vobj{eta}.clnp.il=new_il;
  new_vobj{eta}.clnp.in=new_in;
  new_vobj{eta}.clnp.ip=new_ip;
  new_vobj{eta}.clnp.optflag=new_optflag;
  new_vobj{eta}.clnp.c=new_c;
  new_vobj{eta}.nu=new_nu;
  new_vobj{eta}.cbar = new_vobj{eta}.clnp.c; %for convenience

  %replace existing names by new names so that overwriting is harder
  new_vobj{eta}.clnp_fn=[new_vobj{eta}.clnp_fn '+'];
  new_vobj{eta}.nu_fn=[new_vobj{eta}.nu_fn '+'];
  new_vobj{eta}.q_fn=[new_vobj{eta}.q_fn '+'];

end
