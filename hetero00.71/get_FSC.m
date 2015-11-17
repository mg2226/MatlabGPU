function [magk,fsc]=get_FSC(onevobjA,onevobjB,minmagk,maxmagk,deltamagk)
%function [magk,fsc]=get_FSC(onevobjA,onevobjB,minmagk,maxmagk,deltamagk)

magk=[minmagk:deltamagk:maxmagk].';

if isempty(magk)
  fsc=[];
  return;
end

%Assume that any set of angular basis functions used in the program
%form an orthonormal system on the surface of the sphere.

%Require that the angular basis functions are the same in both models.
%This can be relaxed, but more general cases are not implemented.
if onevobjA.BasisFunctionType~=onevobjB.BasisFunctionType
  error('get_FSC: onevobjA.BasisFunctionType %d ~= onevobjB.BasisFunctionType %d\n', ...
        onevobjA.BasisFunctionType,onevobjB.BasisFunctionType);
end

R1A=onevobjA.R1;
R2A=onevobjA.R2;
llA=onevobjA.clnp.il;
mmA=onevobjA.clnp.in; %Spherical harmonics notation of 'm' not icosahedral harmonics notation of 'n'.
ppA=onevobjA.clnp.ip;
cA=onevobjA.clnp.c;

R1B=onevobjB.R1;
R2B=onevobjB.R2;
llB=onevobjB.clnp.il;
mmB=onevobjB.clnp.in; %Spherical harmonics notation of 'm' not icosahedral harmonics notation of 'n'.
ppB=onevobjB.clnp.ip;
cB=onevobjB.clnp.c;

%Cannot use the Htable that may be a part of either virusobj because
%those values of $|\vec k|$ are unrelated to the values stored in
%magk.  The following code is modified from
%virusobj_set_2Dreciprocal.m .

if R1A>=0.0
  error('get_FSC: Class A, only R1<0.0 (not R1>=0.0) is implemented. R1 %g\n',R1A);
end
[unique_lpA, map_lp2uniqueA, map_unique2lpA] =unique([llA ppA], 'rows');
HtableA=set_H0table_c(unique_lpA(:,1), unique_lpA(:,2), magk, zeros(size(magk)), R2A);

if R1B>=0.0
  error('get_FSC: Class B, only R1<0.0 (not R1>=0.0) is implemented. R1 %g\n',R1B);
end
[unique_lpB, map_lp2uniqueB, map_unique2lpB] =unique([llB ppB], 'rows');
HtableB=set_H0table_c(unique_lpB(:,1), unique_lpB(:,2), magk, zeros(size(magk)), R2B);

SAB=zeros(size(magk));
SAA=zeros(size(magk));
SBB=zeros(size(magk));

%Zhye Yin JSB 2003 Eq. 23.  The l,n,p sums are combined in the ii loop.  The p^\prime sum is the jj loop.

%Compute SAA and SAB
for ii=1:length(llA)

  l=llA(ii);
  m=mmA(ii);
  Hlp=HtableA(:,map_unique2lpA(ii));

  %Find the weights in the B class that have these l,m and any p^\prime
  indices=find( (llB==l) & (mmB==m) );

%  fprintf(1,'get_FSC: ii1 ii %d l %d m %d size(Hlp) %d %d size(indices from B) %d %d\n',ii,l,m,size(Hlp),size(indices));

  if ~isempty(indices)
    for jj=indices'

      Hlpprime=HtableB(:,map_unique2lpB(jj));

%      fprintf(1,'get_FSC: ii1-jj1 ii %d l %d m %d size(Hlp) %d %d jj %d size(Hlpprime) %d %d\n',ii,l,m,size(Hlp),jj,size(Hlpprime));

      SAB=SAB + cA(ii).*cB(jj).*Hlp.*Hlpprime;
    end
  end

  %Find the weights in the A class that have these l,m and any p^\prime
  indices=find( (llA==l) & (mmA==m) );

  if ~isempty(indices)
    for jj=indices'

      Hlpprime=HtableA(:,map_unique2lpA(jj));

%      fprintf(1,'get_FSC: ii1-jj2 ii %d l %d m %d size(Hlp) %d %d jj %d size(Hlpprime) %d %d\n',ii,l,m,size(Hlp),jj,size(Hlpprime));

      SAA=SAA + cA(ii).*cA(jj).*Hlp.*Hlpprime;
    end
  end

end %close for ii=1:length(llA)

%Compute SBB

for ii=1:length(llB)

  l=llB(ii);
  m=mmB(ii);
  Hlp=HtableB(:,map_unique2lpB(ii));

%  fprintf(1,'get_FSC: ii2 ii %d l %d m %d size(Hlp) %d %d\n',ii,l,m,size(Hlp));

  %Find the weights in the B class that have these l,m and any p^\prime
  indices=find( (llB==l) & (mmB==m) );

  if ~isempty(indices)
    for jj=indices'

      Hlpprime=HtableB(:,map_unique2lpB(jj));

%      fprintf(1,'get_FSC: ii2-jj1 ii %d l %d m %d size(Hlp) %d %d jj %d size(Hlpprime) %d %d\n',ii,l,m,size(Hlp),jj,size(Hlpprime));

      SBB=SBB + cB(ii).*cB(jj).*Hlp.*Hlpprime;
    end
  end

end

fsc=[SAB SAA SBB SAB./sqrt(SAA.*SBB)];
