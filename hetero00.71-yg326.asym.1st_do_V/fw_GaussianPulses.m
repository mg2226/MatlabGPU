function [y_nonoise,img,truevalues,pixelnoisevar]=fw_GaussianPulses(SNR,GaussianPulses,deltachi,Nv,NT,Na,Nb,rule)
%function [y_nonoise,img,truevalues,pixelnoisevar]=fw_GaussianPulses(SNR,GaussianPulses,deltachi,Nv,NT,Na,Nb,rule)
%GaussianPulses{1:Neta}.q=class probability
%GaussianPulses{1:Neta}.gainminmax(1:2)=range of uniform pdf on gain for this class which is the sole source of continuous heterogeneity
%GaussianPulses{1:Neta}.means(1:Npulse,1:3)=mean vector for each pulse in the class
%GaussianPulses{1:Neta}.covars(1:Npulse,1:3,1:3)=covariance matrix for each pulse in the class

%Can do only one image per virus particle, no tilt series.
if NT~=1
  error('fw_GaussianPulses: NT %d\n',NT);
end

WEIGHT_INDEX=6; %which column of rule(:,:) has the weights

Neta=length(GaussianPulses);

%Noise-free 2-D reciprocal-space images in the format expected by the estimation functions.
%y_nonoise=zeros(2*size(vkminimalset,1),Nv);
y_nonoise=[];
fprintf(1,'fw_GaussianPulses: y_nonoise is not implemented.  Could be done by fft or by symbolic formula based on characteristic function of a Gaussian pdf.\n')

%Noisy 2-D real-space images.
img=cell(Nv,1);

%Prepare for random selection of class.
q_all=zeros(Neta,1);
for eta=1:Neta
  q_all(eta)=GaussianPulses{eta}.q;
end
cumq=cumsum(q_all);
assert(cumq(end)==1); %the sum of all class probabilities should equal to 1.

%Prepare for random selection of projection direction and origin offset.
%While it would be possible to use a class-specific rule, that is not currently implemented.
cumtheta=cumsum(rule(:,WEIGHT_INDEX),1)/sum(rule(:,WEIGHT_INDEX),1);

%Prepare for evaluating the 2-D Gaussian for the real-space image.
xa=[-floor(Na/2):-floor(Na/2)+Na-1]*deltachi(1);
xb=[-floor(Nb/2):-floor(Nb/2)+Nb-1]*deltachi(2);
[Xa,Xb]=ndgrid(xa,xb);
vecXa=Xa(:);
vecXb=Xb(:);
allpixels=[vecXa' ; vecXb']; %should be 2 x Na*Nb
assert(size(allpixels,1)==2,'fw_GaussianPulses: allpixels: first dimension');
assert(size(allpixels,2)==Na*Nb,'fw_GaussianPulses: allpixels: second dimension');

truevalues=zeros(Nv,2); %truevalues(:,1)=true eta, truevalues(:,2)=true index into rule
mean_all_pixels_inrealspace=0.0;
for ii=1:Nv
  %Determine the realization of the class.
  eta=find((cumq>=rand(1)), 1);
  truevalues(ii,1)=eta;
  %Determine the realization of the projection direction and origin offset.
  itheta=find((cumtheta>=rand(1)), 1);
  truevalues(ii,2)=itheta;
  Rabc=euler2R(rule(itheta,1:3));
  %Determine the gain
  gain=GaussianPulses{eta}.gainminmax(1)+(GaussianPulses{eta}.gainminmax(2)-GaussianPulses{eta}.gainminmax(1))*rand();
  %Sum the contributions of each pulse
  totalimage=zeros(Na,Nb);
  for np=1:size(GaussianPulses{eta}.means,1)
    %Compute the rotated Gaussian
    rmean=Rabc*(GaussianPulses{eta}.means(np,:)');
    rcovar=Rabc*squeeze(GaussianPulses{eta}.covars(np,:,:))*(Rabc.');
    %Project in the z direction = get marginal on x-y = extract subblock
    mean2=rmean(1:2);
    covar2=rcovar(1:2,1:2);
    %Evaluate this Gaussian at all pixels
    normalizer=1/(2*pi*sqrt(det(covar2)));
    tmp=bsxfun(@minus,allpixels,mean2);
    totalimage=totalimage+reshape(normalizer.*exp(-0.5.*sum(tmp.*(covar2\tmp))),Na,Nb);
  end
  img{ii}=totalimage;
  mean_all_pixels_inrealspace=mean_all_pixels_inrealspace+sum(sum(img{ii}));
end
mean_all_pixels_inrealspace=mean_all_pixels_inrealspace/(Nv*Na*Nb);
fprintf(1,'fw_GaussianPulses: mean_all_pixels_inrealspace %g\n',mean_all_pixels_inrealspace);

if SNR==Inf
  fprintf(1,'fw_GaussianPulses: SNR==Inf, will return noise-free 2-D real-space images\n');
  return;
end

%Compute covariance
cov_all_pixels_inrealspace=0.0;
for ii=1:Nv
  cov_all_pixels_inrealspace=cov_all_pixels_inrealspace+sum(sum((img{ii}-mean_all_pixels_inrealspace).^2));
end
cov_all_pixels_inrealspace=cov_all_pixels_inrealspace/(Nv*Na*Nb-1);
fprintf(1,'fw_GaussianPulses: cov_all_pixels_inrealspace %g\n',cov_all_pixels_inrealspace);
fprintf(1,'fw_GaussianPulses: sqrt(cov_all_pixels_inrealspace) %g\n',sqrt(cov_all_pixels_inrealspace));

%Compute the maximum
max_all_pixels_inrealspace=0.0;
for ii=1:Nv
  tmp=max(max(img{ii}));
  max_all_pixels_inrealspace=max(max_all_pixels_inrealspace,tmp);
end
fprintf(1,'fw_GaussianPulses: max_all_pixels_inrealspace %g\n',max_all_pixels_inrealspace);

%Compute the standard deviation of the noise to be added.
%The is a SNR on standard deviation not on variance.
stddev_of_noise=sqrt(cov_all_pixels_inrealspace)/SNR;
fprintf(1,'fw_GaussianPulses: cov_all_pixels_inrealspace %g SNR %g stddev_of_noise %g\n',cov_all_pixels_inrealspace,SNR,stddev_of_noise);
pixelnoisevar=stddev_of_noise.^2;
fprintf(1,'fw_GaussianPulses: pixelnoisevar %g\n',pixelnoisevar);

%Add noise.
for ii=1:Nv
  img{ii}=img{ii} + stddev_of_noise.*randn(size(img{ii}));
end
