function [q,count]=rd_q(fn)
%function [q,count]=rd_q(fn)

fid=fopen(fn);
if fid==-1
  error('rd_q: fopen returned fid -1\n');
end

[q,count]=fscanf(fid,'%g');
if count~=1
  error('rd_q: fscanf returned count %d ~= 1\n',count);
end

if q<0
  error('rd_q: q %g < 0\n',q);
end
