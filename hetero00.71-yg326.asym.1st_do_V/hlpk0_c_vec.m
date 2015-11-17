function y=hlpk0_c_vec(l, p, k, root, rmax)
%function y=hlpk0_c_vec(l, p, k, root, rmax)
%
% Translated from hlpk0.c
%
% This function computes the same results as the C version of hlpk0.c but
% the correct version should flip the sign of the "empirical" variable if
% fabs(df) == -df.
% See
% /home/hours/yz247/Research/FromSeunghee/EM.helix/MLE.serial.motif.known/s
% rc/hlpk0.c
%
% Y. Zheng
% 4/03/2008
%

EPS=1.0e-6;
normhlpk0 = sqrt(2) * sqrt(rmax) * rmax;
tmp = (2.0 * pi .* k) * rmax;
empirical = 4.0 * pi;
[sj, sjp_junk, sy_junk, syp_junk] = sphbes_vec(l, tmp);
numer = normhlpk0 * root(l+1,p) .* sj;
denom = tmp .* tmp - root(l+1,p)^2;

ind = find (abs(numer)<=EPS & abs(denom)<=EPS);
if ~isempty(ind)
  % /* use L'Hospital's Rule */
  numer(ind) = normhlpk0 * root(l+1,p) * sjp_junk(ind);
  denom(ind) = 2.0 * tmp(ind);
  fprintf(2,'hlpk0_c_vec: l %d p %d tmp-root %g, numer %g denom %g by L''Hospital''s Rule\n',...
    l, p, tmp - root(l+1,p), numer, denom);
  
  if ~isempty(find(abs(numer)<=EPS & abs(denom)<=EPS,1))
    %/* should never reach this statement */
    fprintf(2, 'hlpk0_c_vec: bad, l %d p %d k %g\n',l,p,k);
    fprintf(2, 'hlpk0_c_vec: bad, EPS %g normhlpk0 %g root(l+1,p) %g sjp_junk %g tmp %g\n',...
      EPS,normhlpk0,root(l+1,p),sjp_junk,tmp);
    error('should never reach this statement!\n');
  end
end
  
y = empirical * numer ./ denom;
