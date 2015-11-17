function plm = plgndr(l, m, x)
%function plm = plgndr(l, m, x)
% 
% Translated from plgndr.c
% See Numerical Recipes Ch. 6.7
%
% x can be a scalar or a vector
%
% Yili Zheng
% 6/25/2008
%

if (m < 0 || m > l || nnz(abs(x) > 1.0))
  error('Bad arguments in routine plgndr\n');
end

pmm = ones(size(x));

if (m > 0) 
  somx2 = sqrt((1.0 - x) .* (1.0 + x));
  fact = 1.0;
  for i= 1:m
    pmm = pmm .* (-fact * somx2);
    fact = fact + 2.0;
  end
end

if (l == m)
  plm = pmm;
  return;
else 
  pmmp1 = x .* (2 * m + 1) .* pmm;
  if (l == (m + 1))
    plm = pmmp1;
    return;
  else
    for ll = m + 2:l
      pll = (x .* (2 * ll - 1) .* pmmp1 - (ll + m - 1) .* pmm) / (ll - m);
      pmm = pmmp1;
      pmmp1 = pll;
    end
    plm = pll;
  end
end

