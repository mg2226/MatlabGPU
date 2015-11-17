function [vk,y]=subset_reciprocalspace(kmax,vkmag,vkminimalset,Imagestack,iiminimalset)
%function [vk,y]=subset_reciprocalspace(kmax,vkmag,vkminimalset,Imagestack,iiminimalset)

%Construct the reciprocal space image data structure for the range of reciprocal space $\|\vec\kappa\|<kmax$ that will be used in this step.

%I believe we continue to separately use both $\Re Y(\vec\kappa=0)$ and $\Im Y(\vec\kappa=0)$ even though $\Im Y(\vec\kappa=0)=0$.

if kmax==-1 %use everything
  vkindex=[1:length(vkmag)];
  vk=vkminimalset;
elseif kmax==-2 %use everything except |\vec\kappa|=0
  vkindex=find(0<vkmag);
  vk=vkminimalset(vkindex,:);
elseif 0<=kmax %general case
  vkindex=find(vkmag<kmax);
  vk=vkminimalset(vkindex,:);
else %error case
  error('subset_reciprocalspace: kmax %g\n',kmax);
end
Ny=size(vk,1);
%Variable y will store each image in consecutive memory; matlab uses column-major storage.
y=zeros(Ny*2,length(Imagestack));
for ii=1:length(Imagestack)
  y(1:2:end,ii)=real(Imagestack{ii}(iiminimalset(vkindex)));
  y(2:2:end,ii)=imag(Imagestack{ii}(iiminimalset(vkindex)));
end
