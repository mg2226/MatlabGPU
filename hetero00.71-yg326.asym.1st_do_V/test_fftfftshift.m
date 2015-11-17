function test_fftfftshift()

N1=10;
N2=11;
test_fftfftshift2(N1,N2);

N1=11;
N2=10;
test_fftfftshift2(N1,N2);

N1=9;
N2=11;
test_fftfftshift2(N1,N2);

N1=10;
N2=12;
test_fftfftshift2(N1,N2);

function test_fftfftshift2(N1,N2)
%function test_fftfftshift2(N1,N2)

fprintf(1,'>>>>>>>>>>N1 %d N2 %d<<<<<<<<<<\n',N1,N2);

imagestack=rand(N1,N2)+sqrt(-1)*rand(N1,N2);
Imagestack=fft2(ifftshift(imagestack));
imagestack2=fftshift(ifft2(Imagestack));

deltastack=imagestack-imagestack2;
minabsdelta=max(abs(deltastack(:)));
maxabsdelta=min(abs(deltastack(:)));
error=sum(abs(deltastack(:))./abs(imagestack(:)))/(N1*N2);
fprintf(1,'minabsdelta %g maxabsdelta %g error %g\n',minabsdelta,maxabsdelta,error);
