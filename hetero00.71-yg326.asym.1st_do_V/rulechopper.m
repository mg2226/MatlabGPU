function newrule=rulechopper()
%function newrule=rulechopper()

oldrule=rd_rule('rule_small_3');
newrule=rulechopper1(oldrule);
fid=fopen('rule_small_3_rulechopper','w');
if fid==-1
  error('rulechopper: fopen\n');
end
fprintf(fid,'%g %g %g %g %g %g\n',newrule');
fclose(fid);

function newrule=rulechopper1(oldrule)
%function newrule=rulechopper1(oldrule)

rng(2920723); %set the pseudo random number seed

Nz=size(oldrule,1);

probabilitykeep=.01;

keepthisabscissa=( rand(Nz,1)<=probabilitykeep );

newrule=oldrule(keepthisabscissa,:);

newrule=newrule.*size(oldrule,1)./size(newrule,1);
