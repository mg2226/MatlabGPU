AsymObj_ic_lmax2pmax5_clnp_eta1
AsymObj_ic_lmax2pmax5_clnp_eta2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmax=0;
pmax=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn=['AsymObj.ic.lmax' num2str(lmax) 'pmax' num2str(pmax) '.q.equalclassprobs.txt'];
fid=fopen(fn,'w');
if fid==-1
  error('cannot open %s\n',fn);
end

fprintf(fid,'0.5\n');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn=['AsymObj.ic.lmax' num2str(lmax) 'pmax' num2str(pmax) '.q.equalclassprobswith4classes.txt'];
fid=fopen(fn,'w');
if fid==-1
  error('cannot open %s\n',fn);
end

fprintf(fid,'0.25\n');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn=['AsymObj.ic.lmax' num2str(lmax) 'pmax' num2str(pmax) '.q.is1.txt'];
fid=fopen(fn,'w');
if fid==-1
  error('cannot open %s\n',fn);
end

fprintf(fid,'1.0\n');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[status,result]=unix('touch AsymObj.ic.lmax0pmax1.nu.homogeneous.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
