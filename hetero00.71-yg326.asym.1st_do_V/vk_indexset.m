function [Ny4minimalset,iiminimalset,vkminimalset,ixminimalset,iyminimalset]=vk_indexset(Na,Nb,deltachia,deltachib)
%function [Ny4minimalset,iiminimalset,vkminimalset,ixminimalset,iyminimalset]=vk_indexset(Na,Nb,deltachia,deltachib)

[kx,ky]=set_k2D(Na,Nb,deltachia,deltachib);
[ix,iy]=set_indices2D(Na,Nb);

minimalset=get_minimal_set(Na,Nb);

[iiminimalset,vkminimalset,ixminimalset,iyminimalset]=extract(minimalset,kx,ky,ix,iy);

Ny4minimalset=2*max(size(find(minimalset==1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iiset,vkset,ixset,iyset]=extract(set,kx,ky,ix,iy)
%function [iiset,vkset,ixset,iyset]=extract(set,kx,ky,ix,iy)

iiset=find(set);

vkset=zeros([size(iiset,1),2]);
vkset(:,1)=kx(iiset);
vkset(:,2)=ky(iiset);

ixset=ix(iiset);
iyset=iy(iiset);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ka,kb]=set_k2D(Na,Nb,deltachia,deltachib)
%function [ka,kb]=set_k2D(Na,Nb,deltachia,deltachib)
%deltachia and deltachib are the sampling intervals in the
%real space image.  Using values of 1 gives pixel-size units.

k1D_a=set_k1D(Na,deltachia);
k1D_b=set_k1D(Nb,deltachib);

ka=zeros(Na,Nb);
for i=1:Nb
  ka(:,i)=k1D_a;
end

kb=zeros(Na,Nb);
for i=1:Na
  kb(i,:)=k1D_b';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k=set_k1D(N,deltachi)
%function k=set_k1D(N,deltachi)

k=zeros(N,1);
for n=1:N %index into Matlab matrix
  k(n)=index2k(n,N,deltachi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ia,ib]=set_indices2D(Na,Nb)
%function [ia,ib]=set_indices2D(Na,Nb)
%modeled on
%function [ka,kb]=set_k2D(Na,Nb,deltachia,deltachib)

v=[1:Na]';
ia=zeros(Na,Nb);
for i=1:Nb
  ia(:,i)=v;
end

v=1:Nb;
ib=zeros(Na,Nb);
for i=1:Na
  ib(i,:)=v;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function minimal_set=get_minimal_set(ndim1,ndim2)
%function minimal_set=get_minimal_set(ndim1,ndim2)
%(ndim1,ndim2)=dimension of rectangular matrix
%the DC sample is at (1,1)

minimal_set=zeros(ndim1,ndim2);
got_it=zeros(ndim1,ndim2);
for x1=1:ndim1 %index into Matlab matrix
  x2=conj_sym_index_1D(x1,ndim1);
  for y1=1:ndim2 %index into Matlab matrix
    y2=conj_sym_index_1D(y1,ndim2);
    if got_it(x2,y2)==0  %have I already got the conjugate symmetric point?
      %no, I don't have it
      minimal_set(x1,y1)=1;
      got_it(x1,y1)=1;
      got_it(x2,y2)=1;
    end
  end
end
min_got_it=min(got_it(:));
max_got_it=max(got_it(:));
if min_got_it~=1 | max_got_it~=1
  fprintf(1,'get_minimal_set: ndim1 %d ndim2 %d min(got_it(:)) %g max(got_it(:)) %g\n', ...
          ndim1,ndim2,min_got_it,max_got_it);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k=index2k(n,N,deltachi)
%function k=index2k(n,N,deltachi)
%n in 1:N is an index into a Matlab vector/matrix

N2=N/2;

if n<1 | n>N
  fprintf(1,'index2k: old n %d N %d\n',n,N);
  n=mod(n-1,N)+1;
  fprintf(1,'index2k: new n %d N %d\n',n,N);
end

nnear0=n-1;
if nnear0>=N2
  nnear0=nnear0-N;
end
          %units of cycles/(units of deltachi) not radians/(units of deltachi)
k=nnear0/(N*deltachi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
