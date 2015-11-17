function expshift=scaling_rule(minexp, maxexp)
%function expshift=scaling_rule(minexp, maxexp)
%
% This is the Matlab version of scaling_rule.c
% Example: expshit = scaling_rule(min(exponent(:), max(exponent(:));
% 
% Yili Zheng, 8/20/2007
%

SAFETYFRAC = 0.9;

if (isscalar(minexp) == 0 || isscalar(maxexp) == 0)
  error('Either minexp or maxexp is not a scalar!\n')
end

% find out the range of exponent that is allowed
minallow = log(realmin);
maxallow = log(realmax);

%{
    Since the program will do exp(+exponent)*mess summed over many
    terms, I am somewhat concerned about overflow even if the exp
    itself is ok.  So I won't use the entire range of possible
    exponents -- I'd rather let more cases than necessary result in
    underflows.  In particular, I will only use
    [minallow,maxallow*SAFETYFRAC].
%}
maxallow = maxallow * SAFETYFRAC;
tmp = (maxallow-minallow) - (maxexp-minexp);
if tmp > 0.0
  % The entire range of exponents will fit, get it centered inthe range
  %   #ifdef PRTINRANGE
  %   fprintf(2, 'inrange %g %g %g %g\n', minallow, maxallow, minexp, maxexp);
  %   #endif
  expshift = minallow+tmp/2.0-minexp;
else
  %{
   The entire range of exponents will not fit, put the largest at
   the top of the range 
  %}
  %fprintf(2, 'outrange %g %g %g %g\n', minallow, maxallow, minexp, maxexp);
  expshift = maxallow - maxexp;
end

