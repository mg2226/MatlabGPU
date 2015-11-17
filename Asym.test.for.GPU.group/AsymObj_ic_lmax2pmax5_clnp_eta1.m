R2=197.4;

rholl=.5;
rhomm=.7;
rhopp=.6;

%set seed
rng(223031);
%lower and upper limits for the uniform generator.
U1=.7;
U2=1.3;

lmax=2;
pmax=5;

vecl=[];
vecm=[];
vecp=[];
vecclnp=[];
for ll=0:lmax
  for mm=0:2*ll
    %the 2nd rand factor is -1 with prob 1/2 and +1 with prob 1/2
    clnp=(rholl^ll)*(rhomm^mm)*(U1 + rand(pmax,1)*(U2-U1)).*(2*((rand(pmax,1)>0.5)-0.5));
    for pp=1:pmax
      vecl=[vecl ; ll];
      vecm=[vecm ; mm];
      vecp=[vecp ; pp];
      vecclnp=[vecclnp ; clnp(pp)*rhopp.^(pp-1)];
    end
  end
end

%scale by 100 just to make the numbers easier to read
vecclnp=100*vecclnp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn=['AsymObj.ic.lmax' num2str(lmax) 'pmax' num2str(pmax) '.clnp.eta1.txt'];
fid=fopen(fn,'w');
if fid==-1
  error('cannot open %s\n',fn);
end

allelements=[vecl vecm vecp vecclnp];

fprintf(fid,'0\n');
fprintf(fid,'-1 %g\n',R2);
fprintf(fid,'%d %d %d 1 %g\n',allelements');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxclnp=max(vecclnp);

stddevFROMmean=.1;
minstddevFROMmaxmean=.01;

vecnu=(stddevFROMmean.*vecclnp).^2;
ii=find(vecnu<(minstddevFROMmaxmean*maxclnp).^2);
vecnu(ii)=(minstddevFROMmaxmean*maxclnp).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn=['AsymObj.ic.lmax' num2str(lmax) 'pmax' num2str(pmax) '.nu.eta1.txt'];
fid=fopen(fn,'w');
if fid==-1
  error('cannot open %s\n',fn);
end

fprintf(fid,'%g\n',vecnu);

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn=['AsymObj.ic.lmax' num2str(lmax) 'pmax' num2str(pmax) '.q.eta1.txt'];
fid=fopen(fn,'w');
if fid==-1
  error('cannot open %s\n',fn);
end

fprintf(fid,'0.6\n');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmax=0;
pmax=1;
fn=['AsymObj.ic.lmax' num2str(lmax) 'pmax' num2str(pmax) '.clnp.eta1.RealSphHarm.c001iszero.txt'];
fid=fopen(fn,'w');
if fid==-1
  error('cannot open %s\n',fn);
end

fprintf(fid,'0\n');
fprintf(fid,'-1 %g\n',R2);
fprintf(fid,'0 0 1 1 0.0\n');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
