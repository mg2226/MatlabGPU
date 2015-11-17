deltatheta=pi/100;
deltaphi=pi/100;
Nxyz=50;
p=1;
j=1;
l=10;
n=-9999; %unused
b_fn='RealBasisFunctionCoeff_RealUniIrrepPhyPiover2_noStructure_1012.txt';
lmax_in_b=45;

[zeroone,samples]=zero_one_cube_from_angularbasisfunc(deltatheta,deltaphi,Nxyz,p,j,l,n,b_fn,lmax_in_b);

minabssamples=min(abs(samples));
indexsmallest=find(samples==minabssamples);
if length(indexsmallest)~=1
  error('samples: length(indexsmallest)~=1');
end

figure;
nz=indexsmallest;
imagesc(zeroone(:,:,nz));
colorbar;
axis equal
title(['zeroone: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);

figure;
nz=indexsmallest + floor(0.3*(size(zeroone,3)-indexsmallest));
imagesc(zeroone(:,:,nz));
colorbar;
axis equal
title(['zeroone: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);

figure;
nz=indexsmallest + floor(0.6*(size(zeroone,3)-indexsmallest));
imagesc(zeroone(:,:,nz));
colorbar;
axis equal
title(['zeroone: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);

figure;
nz=indexsmallest + floor(0.9*(size(zeroone,3)-indexsmallest));
imagesc(zeroone(:,:,nz));
colorbar;
axis equal
title(['zeroone: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);
