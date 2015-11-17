function qNEW=applyclassproblimit(q,MinimumClassProb)
%function qNEW=applyclassproblimit(q,MinimumClassProb)

if length(q)==1 || MinimumClassProb==0.0
  qNEW=q;
  return;
end

toosmallindices=find(q<MinimumClassProb);

if length(toosmallindices)==0
  fprintf(1,'applyclassproblimit: no change\n');
  qNEW=q;
  return;
end

sumoftoosmall=sum(q(toosmallindices));
deltaprob=length(toosmallindices)*MinimumClassProb - sumoftoosmall;
amount2subtract=deltaprob/(length(q)-length(toosmallindices));

fprintf(1,'applyclassproblimit: q:');
fprintf(1,' %g',q);
fprintf(1,'\n');

bigenoughindices=find(q>=MinimumClassProb);

qNEW=zeros(size(q));
qNEW(toosmallindices)=MinimumClassProb;
qNEW(bigenoughindices)=q(bigenoughindices)-amount2subtract;

fprintf(1,'applyclassproblimit: MinimumClassProb %g qNEW:',MinimumClassProb);
fprintf(1,' %g',qNEW);
fprintf(1,'\n');

if length(find(qNEW<0.0))~=0
  error('applyclassproblimit: qNEW has negative values\n');
end

if sum(qNEW)-1.0>2.0*length(q)*eps*100 %2015-07-13 pd83 inserted factor of 100
  fprintf(1,'applyclassproblimit: qNEW not normalized: qNEW\n');
  qNEW
  sum(qNEW)-1.0
  2.0*length(q)*eps
  error('applyclassproblimit: qNEW not normalized\n');
end
