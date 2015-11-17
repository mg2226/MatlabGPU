function test_FFT_pt2(N1,N2,delta1,delta2,freqintegers,vk,y)
%function test_FFT_pt2(N1,N2,delta1,delta2,freqintegers,vk,y)
%N1,N2: dimensions of the image in pixels
%delta1,delta2: sampling interval of the image
%freqintegers[1:number_of_images,1:2]: integers that describe the 2-D freq vector
%vk[Ny,1:2]: 2-D image frequency vectors
%y[2*Ny,1:Nv]: reciprocal space images as vectors, real part then imaginary part

fprintf(1,'test_FFT_pt2: size(vk) %d %d\n',size(vk,1),size(vk,2));
fprintf(1,'test_FFT_pt2: size(y) %d %d\n',size(y,1),size(y,2));

for ii=1:size(freqintegers,1)
  cmplx_y=complex(y(1:2:end,ii),y(2:2:end,ii));
  largest=max(abs(cmplx_y));
  indices=find(abs(cmplx_y)>1.0e-8*largest);
  if length(indices)>1
    fprintf(1,'test_FFT_pt2: length(indices) %d\n',length(indices));
    error('test_FFT_pt2');
  end
  vkappa=[freqintegers(ii,1)/(N1*delta1) freqintegers(ii,2)/(N2*delta2)];
  fprintf(1,'test_FFT_pt2: ii %d largest %g indices %d vkappa %g %g vk %g %g\n',ii,largest,indices,vkappa(1),vkappa(2),vk(indices,1),vk(indices,2));
end
