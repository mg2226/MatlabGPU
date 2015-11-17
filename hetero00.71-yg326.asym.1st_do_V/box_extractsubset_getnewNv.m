function newNv=box_extractsubset_getnewNv(maxNv,a,b,Nv)
%function newNv=box_extractsubset_getnewNv(maxNv,a,b,Nv)

if maxNv<0 | fix(maxNv)~=maxNv
  error('box_main: box_extractsubset_getnewNv: bad maxNv %g\n',maxNv)
end
if a<=0 | fix(a)~=a
  error('box_main: box_extractsubset_getnewNv: bad a %g\n',a)
end
if a>Nv
  error('box_main: box_extractsubset_getnewNv: a %d > Nv %d\n',a,Nv);
end
if b<=0 | fix(b)~=b
  error('box_main: box_extractsubset_getnewNv: bad b %g\n',b)
end

%require a>0, b>0.
%images used will be a + b*(n-1) for n=1:newNv so the first image is
%'a' and the last image is a + b*(newNv-1).
%find the largest integer newNv such that a + b*(newNv-1) <= Nv.
newNv=fix((Nv-a)/b + 1);
newNv=min(newNv,maxNv);
