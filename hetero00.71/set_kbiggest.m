function kbiggest=set_kbiggest(N,deltachi)
%function kbiggest=set_kbiggest(N,deltachi)
%N is the number of samples.
%Matlab would index 1:N.
%The DSP tradition would be to index 0:N-1.
%deltachi is the sampling interval.

if fix(N)~=N
  error('set_kbiggest: N %g fix(N) %d\n',N,fix(N));
end
if N<=0
  error('set_kbiggest: N %d <= 0\n',N);
end
if deltachi<=0.0
  error('set_kbiggest: deltachi %g <= 0\n',deltachi);
end

if fix(N/2)~=N/2
  %N is odd
  nbiggest=(N-1)/2;
else
  %N is even
  nbiggest=N/2;
end

kbiggest=nbiggest/(N*deltachi);
