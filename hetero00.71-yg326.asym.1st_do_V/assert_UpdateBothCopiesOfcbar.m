function assert_UpdateBothCopiesOfcbar(vobj)
%function assert_UpdateBothCopiesOfcbar(vobj)

Neta=length(vobj);

for eta=1:Neta
  if any( (vobj{eta}.cbar-vobj{eta}.clnp.c)~=0.0 )
    error('assert_UpdateBothCopiesOfcbar: eta %d vobj{eta}.cbar ~= vobj{eta}.clnp.c\n',eta);
  end
end
