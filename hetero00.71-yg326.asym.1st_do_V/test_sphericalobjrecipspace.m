function [vecPbysymbolic,ksymbolic,vecPbyFFT,kFFT]=test_sphericalobjrecipspace(rho0,Ro,k1,deltak,k2,N123,deltax123)
%function [vecPbysymbolic,ksymbolic,vecPbyFFT,kFFT]=test_sphericalobjrecipspace(rho0,Ro,k1,deltak,k2,N123,deltax123)

%symbolic
fprintf(1,'test_sphericalobjrecipspace: >>>>>>>>>>symbolic calculations<<<<<<<<<<\n');

flag=0;
if k1==0.0
  flag=1;
  fprintf(1,'test_sphericalobjrecipspace: k1==0.0\n');
  k1=deltak;
  P0=(4/3)*pi*Ro^3*rho0;
  fprintf(1,'test_sphericalobjrecipspace: P0 %g\n',P0);
end
ksymbolic=k1:deltak:k2;
tmp=2*pi*ksymbolic*Ro;

[sj,sjp_junk,sy_junk,syp_junk]=sphbes_vec(1,tmp);
vecPbysymbolic=(rho0*4*pi*Ro^3).*sj./tmp;

if flag==1
  ksymbolic=[0 ksymbolic];
  vecPbysymbolic=[P0 vecPbysymbolic];
end

fprintf(1,'test_sphericalobjrecipspace: symbolic %g %g ; %g %g\n', ...
        ksymbolic(1),vecPbysymbolic(1),ksymbolic(2),vecPbysymbolic(2));

%FFT
fprintf(1,'test_sphericalobjrecipspace: >>>>>>>>>>FFT calculations<<<<<<<<<<\n');

range123=[-N123:N123];
[x,y,z]=ndgrid(range123,range123,range123);
rsqu=x.^2+y.^2+z.^2;
ii=find(rsqu<=(Ro/deltax123)^2);
rho=zeros(size(x));
rho(ii)=rho0;

fprintf(1,'test_sphericalobjrecipspace: size(range123) %d %d\n',size(range123));
fprintf(1,'test_sphericalobjrecipspace: size(rho) %d %d %d\n',size(rho));
fprintf(1,'test_sphericalobjrecipspace: max(rho(:)) %g\n',max(rho(:)));
fprintf(1,'test_sphericalobjrecipspace: min(rho(:)) %g\n',min(rho(:)));
sumofrhocubetimesvoxel=(deltax123^3)*sum(rho(:));
fprintf(1,'test_sphericalobjrecipspace: sumofrhocubetimesvoxel %g\n',sumofrhocubetimesvoxel);
figure;
imagesc(rho(:,:,floor(end/2)));
axis equal
title('rho(:,:,floor(end/2))');

PbyFFT=fftn(ifftshift(rho)); %as done in hetero for command 'realrecip_2DFFT'
vecPbyFFT=(deltax123^3).*PbyFFT(1:N123,1,1);
kFFT=[0:N123-1]'./((2*N123+1).*deltax123); %frequencies along the first of three axes

fprintf(1,'test_sphericalobjrecipspace: size(vecPbyFFT) %d %d\n',size(vecPbyFFT));

maxabsPbyFFT=max(abs(PbyFFT(:)));
fprintf(1,'test_sphericalobjrecipspace: maxabsPbyFFT %g\n',maxabsPbyFFT);
[i1,i2,i3]=ind2sub(size(PbyFFT),find(abs(PbyFFT)==maxabsPbyFFT));
fprintf(1,'test_sphericalobjrecipspace: the indices of maxabsPbyFFT %d %d %d\n',i1,i2,i3);
maxabsPbyFFTtimesvoxel=(deltax123^3)*maxabsPbyFFT;
fprintf(1,'test_sphericalobjrecipspace: maxabsPbyFFTtimesvoxel %g\n',maxabsPbyFFTtimesvoxel);
fprintf(1,'test_sphericalobjrecipspace: vecPbyFFT(1) %g + i %g\n',real(vecPbyFFT(1)),imag(vecPbyFFT(1)));

%FFT on a constant object
fprintf(1,'test_sphericalobjrecipspace: >>>>>>>>>>FFT on a constant object<<<<<<<<<<\n');

%summary of the following results:
%The DC term in the DFT is located at (1,1,1).
%The DC value is the sum of all the samples.
N1=10;
N2=11;
N3=12;
x=ones(N1,N2,N3);
X=fftn(ifftshift(x));
maxabsX=max(abs(X(:)));
[i1,i2,i3]=ind2sub(size(X),find(abs(X)==maxabsX));
fprintf(1,'test_sphericalobjrecipspace: max(abs(X)) is X(%d,%d,%d)=%g\n',i1,i2,i3,X(i1,i2,i3));
fprintf(1,'test_sphericalobjrecipspace: N1*N2*N3 %d\n',N1*N2*N3);
