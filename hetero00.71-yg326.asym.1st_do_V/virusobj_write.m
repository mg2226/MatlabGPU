function virusobj_write(fn_clnp,fn_nu,fn_q,vobj)
%function virusobj_write(fn_clnp,fn_nu,fn_q,vobj)

for eta=1:length(vobj)

  fid=fopen(fn_clnp{eta},'w');
  if fid==-1
    error('virusobj_write: eta %d fopen fn_clnp %s\n',eta,fn_clnp{eta});
  end
  fprintf(fid,'%d\n',vobj{eta}.BasisFunctionType);
  fprintf(fid,'%g %g\n',vobj{eta}.R1,vobj{eta}.R2);
  c=[vobj{eta}.clnp.il vobj{eta}.clnp.in vobj{eta}.clnp.ip vobj{eta}.clnp.optflag vobj{eta}.clnp.c];
  fprintf(fid,'%d %d %d %d %g\n',c.');
  status=fclose(fid);
  if status~=0
    error('virusobj_write: eta %d fclose fn_clnp %s\n',eta,fn_clnp{eta});
  end

  if ~isempty(vobj{eta}.nu)

    fid=fopen(fn_nu{eta},'w');
    if fid==-1
      error('virusobj_write: eta %d fopen fn_nu %s\n',eta,fn_nu{eta});
    end
    fprintf(fid,'%g\n',vobj{eta}.nu.');
    status=fclose(fid);
    if status~=0
      error('virusobj_write: eta %d fclose fn_nu %s\n',eta,fn_nu{eta});
    end

  else

    unix(['touch ' fn_nu{eta}]);

  end

  fid=fopen(fn_q{eta},'w');
  if fid==-1
    error('virusobj_write: eta %d fopen fn_q %s\n',eta,fn_q{eta});
  end
  fprintf(fid,'%g\n',vobj{eta}.q);
  status=fclose(fid);
  if status~=0
    error('virusobj_write: eta %d fclose fn_q %s\n',eta,fn_q{eta});
  end

end
