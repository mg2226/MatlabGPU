function rule=quadrule4asym1(a0,a1,Na,b0,b1,Nb,c0,c1,Nc)
%function rule=quadrule4asym1(a0,a1,Na,b0,b1,Nb,c0,c1,Nc)
%quadrature rule 1 for Euler angles for asymmetric structures

assert(a0>=0,'quadrule4asym1: a0 negative');
assert(a1>=a0,'quadrule4asym1: a0,a1 order');
assert(a1<=2*pi,'quadrule4asym1: a1>2*pi');

assert(b0>=0,'quadrule4asym1: b0 negative');
assert(b1>=b0,'quadrule4asym1: b0,b1 order');
assert(b1<=2*pi,'quadrule4asym1: b1>2*pi');

assert(c0>=0,'quadrule4asym1: c0 negative');
assert(c1>=c0,'quadrule4asym1: c0,c1 order');
assert(c1<=2*pi,'quadrule4asym1: c1>2*pi');

%uniform integration rule for "d\phi"
astep=(a1-a0)/Na;
aabscissa=a0+astep*[0:Na-1]';
aweights=astep*ones(Na,1);

%Gauss-Legendre integration rule for "d\theta"
[babscissa,bweights]=lgwt(Nb,b0,b1);

%To get d\Omega (Doerschuk and Johnson, IEEE T-IT, 2000; Section II) from a
%and b, multiply total weight aweights(na)*bweights(nb) by sin(babscissa).
%Easiest way to do this is to scale the bweights by sin(babscissa)

%uniform integration rule for "d\psi"
cstep=(c1-c0)/Nc;
cabscissa=c0+cstep*[0:Nc-1]';
cweights=cstep*ones(Nc,1);

%To get p_S (Doerschuk and Johnson, IEEE T-IT, 2000; Section II) from \Omega
%(from a and b), multiply weights by 1/(4*pi)
%To get p_C (Doerschuk and Johnson, IEEE T-IT, 2000; Section II) from c,
%multiply weights by 1/(2*pi)
%Easiest way to do both is to scale the bweights by 1/(8*pi^2).

%Scale the bweights for everything:
bweights=bweights.*sin(babscissa)/(8*pi^2);

%alpha varies slowest, beta is intermediate, and gamma varies fastest
rule=zeros(Na*Nb*Nc,6);
nabc=1;
for na=1:Na
  for nb=1:Nb
    for nc=1:Nc
      rule(nabc,1)=aabscissa(na);
      rule(nabc,2)=babscissa(nb);
      rule(nabc,3)=cabscissa(nc);
      %rule(nabc,4)=saved for first origin offset translation;
      %rule(nabc,5)=saved for first origin offset translation;
      rule(nabc,6)=aweights(na)*bweights(nb)*cweights(nc);
      nabc=nabc+1;
    end
  end
end
