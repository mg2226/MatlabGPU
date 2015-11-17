function imagestack=test_FFT_pt1(N1,N2,delta1,delta2,freqintegers)
%function imagestack=test_FFT_pt1(N1,N2,delta1,delta2,freqintegers)
%N1,N2: dimensions of the image in pixels
%delta1,delta2: sampling interval of the image
%freqintegers[1:number_of_images,1:2]: integers that describe the 2-D freq vector

vecx1=delta1.*[1:N1];
vecx2=delta2.*[1:N2];
[x1,x2]=ndgrid(vecx1,vecx2);

imagestack=cell(size(freqintegers,1),1);

for ii=1:size(freqintegers,1)
  vkappa=[freqintegers(ii,1)/(N1*delta1) freqintegers(ii,2)/(N2*delta2)];
  imagestack{ii}=cos(2*pi*(vkappa(1).*x1+vkappa(2).*x2));
end
