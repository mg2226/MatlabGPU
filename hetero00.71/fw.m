function [y_nonoise,img,truevalues,pixelnoisevar]=fw(SNR,vobj,vkminimalset,Nv,NT,Na,Nb,rule,all_tilde_b,ixminimalset,iyminimalset)
%function [y_nonoise,img,truevalues,pixelnoisevar]=fw(SNR,vobj,vkminimalset,Nv,NT,Na,Nb,rule,all_tilde_b,ixminimalset,iyminimalset)

%Can do only one image per virus particle, no tilt series.
if NT~=1
  error('fw: NT %d\n',NT);
end

if ndims(vkminimalset)~=2 || size(vkminimalset,2)~=2
  error('fw: dimensions of vkminimalset\n');
end

WEIGHT_INDEX=6; %which column of rule(:,:) has the weights

Neta=length(vobj);

%Noise-free 2-D reciprocal-space images in the format expected by the estimation functions.
y_nonoise=zeros(2*size(vkminimalset,1),Nv);
%Noisy 2-D real-space images.
img=cell(Nv,1);

%Prepare for random selection of class.
q_all=zeros(Neta,1);
for eta=1:Neta
  q_all(eta)=vobj{eta}.q;
end
cumq=cumsum(q_all);
assert(cumq(end)==1); %the sum of all class probabilities should equal to 1.

%Prepare for random selection of projection direction and origin offset.
%While it would be possible to use a class-specific rule, that is not currently implemented.
cumtheta=cumsum(rule(:,WEIGHT_INDEX),1)/sum(rule(:,WEIGHT_INDEX),1);

%Compute the Cholesky factorization of V (the weight covariance matrix).
choleskyfactor=cell(Neta,1);
for eta=1:Neta
  if isempty(vobj{eta}.nu)
    choleskyfactor{eta}=[];
  else
    choleskyfactor{eta}=sqrt(vobj{eta}.nu);
  end
end

%Compute 2-D reciprocal-space images by y=Lc.
truevalues=zeros(Nv,2); %truevalues(:,1)=true eta, truevalues(:,2)=true index into rule
mean_all_pixels=0.0;

for ii=1:Nv
  %Determine the realization of the class.
  eta=find((cumq>=rand(1)), 1);
  
  truevalues(ii,1)=eta;
  %Determine the realization of the projection direction and origin offset.
  itheta=find((cumtheta>=rand(1)), 1);
  truevalues(ii,2)=itheta;
  %Determine the realization of the coefficients.
  Nc=length(vobj{eta}.clnp.c);
  c=vobj{eta}.clnp.c;
  if ~isempty(vobj{eta}.nu)
    c=c+choleskyfactor{eta}.*randn(Nc,1);
  end
  %Compute L
  Rabc=euler2R(rule(itheta,1:3));
  L=setL_nord(Rabc, vobj{eta}.clnp.il, vobj{eta}.clnp.in, vkminimalset, vobj{eta}.Htable, vobj{eta}.map_unique2lp, all_tilde_b{eta}.tilde_b);
  if ~isreal(L)
    error('fw: ii %d L is not real\n',ii);
  end

  if (rank(L) < Nc)
    fprintf(2,'L is not full rank! ii %d rank(L) %d Nc %d.\n',ii,rank(L),Nc)
  end
  %Compute 2-D reciprocal-space images in the format expected by the estimation functions.
  y_nonoise(:,ii)=L*c;
  %Compute 2-D real-space images as 2-D arrays.
  ImgAsRealVector=y_nonoise(:,ii);
  ImgAsComplexVector=(ImgAsRealVector(1:2:end)+sqrt(-1).*ImgAsRealVector(2:2:end));
  ImgAs2DComplexImage=recip_Cvec2Cmat_conj_sym(ImgAsComplexVector,Na,Nb,ixminimalset,iyminimalset);
  img{ii}=fftshift(ifft2(ImgAs2DComplexImage));
  mean_all_pixels=sum(sum(img{ii}));
end
mean_all_pixels=mean_all_pixels/(Nv*Na*Nb);

if SNR==Inf
  fprintf(1,'fw: SNR==Inf, will return noise-free 2-D real-space images\n');
  return;
end

%Compute covariance
cov_all_pixels=0.0;
for ii=1:Nv
  cov_all_pixels=cov_all_pixels+sum(sum((img{ii}-mean_all_pixels).^2));
end
cov_all_pixels=cov_all_pixels/(Nv*Na*Nb-1);

%Compute the standard deviation of the noise to be added.
%The is a SNR on standard deviation not on variance.
stddev_of_noise=sqrt(cov_all_pixels)/SNR;
fprintf(1,'fw: cov_all_pixels %g SNR %g stddev_of_noise %g\n',cov_all_pixels,SNR,stddev_of_noise);
pixelnoisevar=stddev_of_noise.^2;
fprintf(1,'fw: pixelnoisevar %g\n',pixelnoisevar);

%Add noise.
for ii=1:Nv
  img{ii}=img{ii} + stddev_of_noise.*randn(size(img{ii}));
end
