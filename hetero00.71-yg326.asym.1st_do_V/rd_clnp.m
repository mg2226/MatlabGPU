function [header1,header2,clnp]=rd_clnp(fn)
%function [header1,header2,clnp]=rd_clnp(fn)

fd=fopen(fn,'r');
if fd==-1
    error('rd_clnp: cannot open fn %s\n',fn);
end

[header1,nparam]=fscanf(fd,'%d\n',1);
if nparam~=1
    fprintf(1,'rd_clnp: bad file header1 nparam %d\n',nparam);
end

[header2,nparam]=fscanf(fd,'%g%g\n',2);
if nparam~=2
    fprintf(1,'rd_clnp: bad file header2 nparam %d\n',nparam);
end

[c,nc]=fscanf(fd,'%g%g%g%d%g\n',inf); %l,n,p,optflag,clnp
if mod(nc,5)~=0
    fprintf(1,'rd_clnp: bad nc %d\n',nc);
end
c=reshape(c,5,nc/5);
c=c';
fclose(fd);

%for real numbers
clnp.il=c(1:1:end,1);
clnp.in=c(1:1:end,2);
clnp.ip=c(1:1:end,3);
clnp.optflag=c(1:1:end,4);
clnp.c = c(1:1:end,5);

%for complex numbers
% clnp.il=c(1:2:end,1);
% clnp.in=c(1:2:end,2);
% clnp.ip=c(1:2:end,3);
% clnp.optflag=c(1:2:end,4);
% clnp.c = c(1:2:end,5)+i*c(2:2:end,5);
