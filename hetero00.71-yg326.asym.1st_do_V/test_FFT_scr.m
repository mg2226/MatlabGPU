%http://www.mathworks.com/matlabcentral/newsreader/view_thread/285244

%http://www.mathworks.com/matlabcentral/fileexchange/25473-why-use-fftshift-fft-fftshift-x----in-matlab-instead-of-fft-x--
%fftshift(fft(ifftshift(sig))) or fftshift(ifft(ifftshift(spectrum))).

N1=20;
N2=21;
delta1=2.74;
delta2=3.24;
vecx1=delta1.*[1:N1];
vecx2=delta2.*[1:N2];
[x1,x2]=ndgrid(vecx1,vecx2);

vkappa=[1/(N1*delta1) 1/(N2*delta2)];

sigma_cos=cos(2*pi*(vkappa(1).*x1+vkappa(2).*x2));
sigma_exp=exp(sqrt(-1)*2*pi*(vkappa(1).*x1+vkappa(2).*x2));

figure;
imagesc(sigma_cos);
axis equal
title('sigma for cos');
colorbar

%abs(sigma_exp)=1 so plot phase
figure;
imagesc(angle(sigma_exp));
axis equal
title('angle(sigma) for exp, |sigma|=1');
colorbar

Sigma_cos=fft2(sigma_cos);

figure;
imagesc(abs(Sigma_cos));
axis equal
title('|Sigma| for cos');
colorbar

absSigma_cos=abs(Sigma_cos);
maxabsSigma_cos=max(absSigma_cos(:));
[ii,jj]=find(absSigma_cos==maxabsSigma_cos)

Sigma_exp=fft2(sigma_exp);

figure;
imagesc(abs(Sigma_exp));
axis equal
title('|Sigma| for exp');
colorbar

absSigma_exp=abs(Sigma_exp);
maxabsSigma_exp=max(absSigma_exp(:));
[ii,jj]=find(absSigma_exp==maxabsSigma_exp)
