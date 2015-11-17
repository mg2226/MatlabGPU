function FHV_create_example()
%function FHV_create_example()

[header1 header2 clnp]=rd_clnp('FHV.lmax10pmax5.clnp.txt');
stddev=mean2stddev(clnp.c,.0025,.02);
fid=fopen('FHV.lmax10pmax5.nu.txt','w');
if fid==-1
  error('FHV_create_example: fopen #1 failed\n');
end
fprintf(fid,'%g\n',stddev);
fclose(fid);

[header1 header2 clnp]=rd_clnp('FHV.lmax10pmax5.clnp.perturbed.txt');
stddev=mean2stddev(clnp.c,.0025,.02);
fid=fopen('FHV.lmax10pmax5.nu.perturbed.txt','w');
if fid==-1
  error('FHV_create_example: fopen #2 failed\n');
end
fprintf(fid,'%g\n',stddev);
fclose(fid);

fid=fopen('FHV.lmax10pmax5.q.txt','w');
if fid==-1
  error('FHV_create_example: fopen #3 failed\n');
end
fprintf(fid,'%g\n',0.6);
fclose(fid);

fid=fopen('FHV.lmax10pmax5.q.perturbed.txt','w');
if fid==-1
  error('FHV_create_example: fopen #4 failed\n');
end
fprintf(fid,'%g\n',0.4);
fclose(fid);

fid=fopen('FHV.lmax10pmax5.q.is1.txt','w');
if fid==-1
  error('FHV_create_example: fopen #5 failed\n');
end
fprintf(fid,'%g\n',1.0);
fclose(fid);
