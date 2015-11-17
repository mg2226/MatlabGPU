function [x, normp]=jroots_c(l,pmax,xmax)
%function [x, normp]=jroots_c(l,pmax,xmax)
% xmax = rtmax
%
% Translated from jroots.c
%
% Y. Zheng
% 2/13/2008
%
% Todo: The normp values may not be correct.  I was unable to verify it.

JMAX=40;
xacc=1.0e-17;
p = 1;
rtn = 0.0;
pi2 = pi/2.0;

x = zeros(pmax,1);
normp = zeros(pmax,1);

while (p <= pmax && rtn <= xmax)
  x1 = rtn + pi2;
  x2 = x1 + pi2;
  for j = 1:JMAX
    [f1,d2,d1,d3] = sphbes(l,x1); %sphbes(l, x1, &f1, &d1, &d2, &d3);
    [f2,d2,d1,d3] = sphbes(l,x2); %sphbes(l, x2, &f2, &d1, &d2, &d3);
    if (f1 * f2 <= 0.)
      break;
    end
    x1 = x2;
    x2 = x1 + pi2;
  end
  rtn = 0.5 * (x1 + x2);
  for j = 1:JMAX
    [f,df,sy,syp] = sphbes(l,rtn); %sphbes(l, rtn, &f, &sy, &df, &syp);
    dx = f / df;
    rtn = rtn - dx;
    if ((x1 - rtn) * (rtn - x2) < 0.0)
      if(f * f1 < 0.)
        x2 = rtn + dx;
        f2 = f;
      else
        x1 = rtn + dx;
        f1 = f;
      end
      rtn = 0.5 * (x1 + x2);
    end
    if (abs(dx) < xacc)
      break;
    end
  end
  x(p) = rtn;
  normp(p) = df / sqrt(2.0);  % /* 10 March 07: sqrt(z*z)=fabs(z) not z */
  % The following fix doesn't work with the clnp's from the current
  % hlpk0_c.
  %normp(p) = abs(df) / sqrt(2.0);  % /* 10 March 07: sqrt(z*z)=fabs(z) not z */
  p = p + 1;
end