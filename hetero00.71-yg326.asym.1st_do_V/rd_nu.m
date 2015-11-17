function [nu,count]=rd_nu(fn)
%function [nu,count]=rd_nu(fn)

fid=fopen(fn);
if fid==-1
  fprintf(1,'rd_nu: fn %s\n',fn);
  error('rd_nu: fopen returned fid -1\n');
end

[nu,count]=fscanf(fid,'%g');
if count<0
  error('rd_nu: fscanf returned count %d < 0\n',count);
end
if count==0
  fprintf(1,'rd_nu: nu has zero length.  Therefore will solve a homogeneous problem.\n');
end

ii=find(nu<0.0);
if length(ii)>0
  error('rd_nu: %d negative elements in nu\n',length(ii));
end
