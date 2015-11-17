function [zeroone,samples]=zero_one_cube_from_angularbasisfunc(deltatheta,deltaphi,Nxyz,p,j,l,n,b_fn,lmax_in_b)
%function [zeroone,samples]=zero_one_cube_from_angularbasisfunc(deltatheta,deltaphi,Nxyz,p,j,l,n,b_fn,lmax_in_b)
%test with
%b_fn='RealBasisFunctionCoeff_RealUniIrrepPhyPiover2_noStructure_1012.txt';
%lmax_in_b=45;
%Peter C. Doerschuk 2014-12-31

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read \cal D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tildeb,space_check,ll,Ipl]=rd_b_general(b_fn,lmax_in_b);
%get the particular basis function that was requested--use j as the row index and ignore n
calDlp=tildeb{l+1,p};
calDlpj=calDlp(j,:); %vector of \cal D coefficients for -l<=m<=+l
if size(calDlpj,1)~=1 || size(calDlpj,2)~=2*l+1
  error(['size(calDlpj) ' num2str(size(calDlj))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute values of I_{p;j,l,n}(\theta,\phi) on a rectangular grid in theta,phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set up the theta,phi grid
thetavalues=[0:deltatheta:pi];
phivalues=[0:deltaphi:2*pi];

%Will need Y_{l,m}(\theta,\phi) for
%the one fixed value of l,
%all values of m (-l<=m<=+l),
%and all values of (\theta,\phi)

%Matlab legendre function provides only has 0<=m<=+l not -l<=m<=+l
%Matlab legendre function writes P_n^m while PCD writes P_{l,m}:
%Matlab m = PCD m = order
%Matlab n = PCD l = degree
%P: row index is m, column index is element of costhetavalues
costhetavalues=cos(thetavalues);
P=legendre(l,costhetavalues);

%k! = factorialvalues(k+1) for k=0,...,2*l
factorialvalues=factorial(0:2*l);
tmp=(2*l+1)/(4*pi);

Ipjln=zeros(length(phivalues),length(thetavalues));
for nphi=1:length(phivalues)
  Y=zeros(size(P));
  expimphi=exp(i.*[0:l].*phivalues(nphi));
  for m=0:l
%   Y_{l,m}(\theta,\phi) for fixed l,phi varying m,theta
    Y(m+1,:)=sqrt(tmp*factorialvalues(l-m+1)/factorialvalues(l+m+1)) ...
             * P(m+1,:)*expimphi(m+1);
  end
% I_{p;j,l,n}(\theta,\phi) for fixed l,phi varying theta
% add the terms for m^\prime=0,...,l for all values of theta
  Ipjln(nphi,:)=calDlpj(l+1:end)*Y;
% add the terms for m^\prime=-l,...,-1 (note reversed order of calDlpj elements)
  signfromm=1-2*mod(1:l,2);
  Ipjln(nphi,:)=Ipjln(nphi,:)+(signfromm.*calDlpj(l:-1:1))*conj(Y(2:end,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test if I_{p;j,l,n}(\theta,\phi) is real, as is expected.
%if not, fix it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isreal(Ipjln)
  fprintf(1,'Ipjln has an imaginary part, but it may be all zeros\n');
else
  fprintf(1,'Ipjln passed ~isreal(Ipjln) test\n');
end
if any(imag(Ipjln(:)))
  fprintf(1,'Ipjln has an imaginary part and at least one element is not zero\n');
  fprintf(1,'will try to fix Ipjln\n');
  notzero=find(abs(Ipjln)~=0);
  fprintf(1,'size(notzero,1) %d prod(size(Ipjln)) %d\n',size(notzero,1),prod(size(Ipjln)))
  err=1-abs(real(Ipjln(notzero)))./abs(Ipjln(notzero)); %will be in 0<=err<=1
  fprintf(1,'min(err(:)) %d\n',min(err(:)));
  fprintf(1,'max(err(:)) %d\n',max(err(:)));
  baderrs=find(abs(err)>1e-12); %perfect arithmetic would not need abs()
  if length(baderrs)>0
    fprintf(1,'length(baderrs) %d\n',length(baderrs))
    error('large imaginary parts\n');
  else
    fprintf(1,'no bad errors, will fix small errors\n');
  end
  tmp=real(Ipjln);
  Ipjln=tmp;
else
  fprintf(1,'Ipjln passed any(imag(Ipjln(:))) test\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scale I_{p;j,l,n}(\theta,\phi) to [0.5,1.0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minIpjln=min(Ipjln(:));
maxIpjln=max(Ipjln(:));
%make Ipjln range from 0 to 1
Ipjln=(Ipjln-minIpjln)/(maxIpjln-minIpjln);
%make Ipjln range from 0.5 to 1
Ipjln=0.5*(1+Ipjln);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot scaled I_{p;j,l,n}(\theta,\phi) on the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure;
%mesh(thetavalues/(2*pi),phivalues/(2*pi),Ipjln) %thetavalues=columns of Ipjln, phivalues=rwos of Ipjln
%title('Ipjln after scaling');
%xlabel('\theta/(2\pi)'); %thetavalues=columns of Ipjln
%ylabel('\phi/(2\pi)'); %phivalues=rows of Ipjln
%fprintf(1,'about to pause\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup the real-space cube and change from rectangular to spherical coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oversize=1.05;
deltaxyz=2*oversize/Nxyz;
samples=[-oversize:deltaxyz:oversize];
minabssamples=min(abs(samples));
indexsmallest=find(samples==minabssamples);
if length(indexsmallest)~=1
  error('samples: length(indexsmallest)~=1');
end
fprintf(1,'min(abs(samples)) %d\n',min(abs(samples)));
[x,y,z]=meshgrid(samples,samples,samples);
%[azimuth,elevation,r]=cart2sph(X,Y,Z);
[phi,theta,radius]=cart2sph(x,y,z);
phiindex=find(phi(:)<0);
phi(phiindex)=2*pi+phi(phiindex);
if ~isfinite(phi(:))
  error('phi fails isfinite');
end
if any(phi(:)<0)
  phiindex=find(phi(:)<0);
  phi(phiindex)/(2*pi)
  error('phi<0');
end
if any(phi(:)>2*pi)
  phiindex=find(phi(:)>2*pi);
  phi(phiindex)/(2*pi)
  error('phi>2*pi');
end
theta=pi/2-theta;
if ~isfinite(theta(:))
  error('theta fails isfinite');
end
if any(theta(:)<0)
  thetaindex=find(theta(:)<0);
  theta(thetaindex)/(2*pi)
  error('theta<0');
end
if any(theta(:)>pi)
  thetaindex=find(theta(:)>pi);
  theta(thetaindex)/(2*pi)
  error('theta>pi');
end
%maxphi=max(phi(:))/(2*pi)
%minphi=min(phi(:))/(2*pi)
%maxtheta=max(theta(:))/(2*pi)
%mintheta=min(theta(:))/(2*pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%round theta,phi onto the grid where Ipjln and represent as matlab indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

itheta2=round(theta/deltatheta)+1; %+1 because matlab arrays start at index 1
if any(itheta2(:)<1)
  badindex=find(itheta2(:)<1);
  itheta2(badindex)
  error('itheta2<1');
end
if any(itheta2(:)==length(thetavalues)+1)
  badindex=find(itheta2(:)==length(thetavalues)+1)
  itheta2(badindex)=length(thetavalues);
end
if any(itheta2(:)>length(thetavalues))
  badindex=find(itheta2(:)>length(thetavalues))
  itheta2(badindex)
  error('itheta2>length(thetavalues)');
end
iphi2=round(phi/deltaphi)+1; %+1 because matlab arrays start at index 1
if any(iphi2(:)<1)
  badindex=find(iphi2(:)<1);
  iphi2(badindex)
  error('iphi2<1');
end
if any(iphi2(:)==length(phivalues)+1)
  badindex=find(iphi2(:)==length(phivalues)+1)
  iphi2(badindex)=length(phivalues);
end
if any(iphi2(:)>length(phivalues))
  badindex=find(iphi2(:)>length(phivalues))
  iphi2(badindex)
  error('iphi2>length(phivalues)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%look up the values at every site in the cube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test code that fails, seems to give a Kronecker product type of answer:
%nrow=3;
%ncol=5;
%x=rand(nrow,ncol)
%mrow=6;
%mcol=7;
%irow=randi(nrow,[mrow mcol])
%icol=randi(ncol,[mrow mcol])
%y=x(irow,icol)
%test code that works:
%nrow=3;
%ncol=5;
%x=rand(nrow,ncol)
%mrow=6;
%mcol=7;
%irow=randi(nrow,[mrow mcol])
%icol=randi(ncol,[mrow mcol])
%linearindex=sub2ind(size(x),irow(:),icol(:));
%y=x(linearindex);
%y=reshape(y,[mrow mcol])

linearindex=sub2ind(size(Ipjln),iphi2(:),itheta2(:));
Ipjln_3D=Ipjln(linearindex);
Ipjln_3D=reshape(Ipjln_3D,size(radius));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make the object whose radius at theta,phi is related to I_{p;j,l,n}(\theta,\phi)
%zeroone=1 if radius at that voxel is less than I_{p;j,l,n}(\theta,\phi) for the \theta,\phi for that voxel.
%zeroone=0 otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zeroone=(radius<Ipjln_3D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make plots of slices throught Ipjln_3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure;
%mesh(samples,samples,Ipjln_3D(:,:,indexsmallest));
%xlabel('x');
%ylabel('y');
%title(['Ipjln\_3D: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at indexsmallest ' num2str(indexsmallest) ' samples(indexsmallest) ' num2str(samples(indexsmallest))]);
%fprintf(1,'about to pause\n');
%pause;

figure;
imagesc(Ipjln_3D(:,:,indexsmallest));
colorbar
axis equal
xlabel('x');
ylabel('y');
title(['Ipjln\_3D: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at indexsmallest ' num2str(indexsmallest) ' samples(indexsmallest) ' num2str(samples(indexsmallest))]);
%fprintf(1,'about to pause\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make plots of slices through the 3-D radius, theta, and phi arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
imagesc(radius(:,:,indexsmallest));
colorbar
title(['radius: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at indexsmallest ' num2str(indexsmallest) ' samples(indexsmallest) ' num2str(samples(indexsmallest))]);
%fprintf(1,'about to pause\n');
%pause;

figure;
imagesc(theta(:,:,indexsmallest));
colorbar;
title(['theta: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at indexsmallest ' num2str(indexsmallest) ' samples(indexsmallest) ' num2str(samples(indexsmallest))]);
%fprintf(1,'about to pause\n');
%pause;

figure;
imagesc(phi(:,:,indexsmallest));
colorbar;
title(['phi: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at indexsmallest ' num2str(indexsmallest) ' samples(indexsmallest) ' num2str(samples(indexsmallest))]);
%fprintf(1,'about to pause\n');
%pause;

nz=indexsmallest + floor(0.3*(size(zeroone,3)-indexsmallest));

figure;
imagesc(radius(:,:,nz));
colorbar
title(['radius: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);
%fprintf(1,'about to pause\n');
%pause;

figure;
imagesc(theta(:,:,nz));
colorbar;
title(['theta: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);
%fprintf(1,'about to pause\n');
%pause;

figure;
imagesc(phi(:,:,nz));
colorbar;
title(['phi: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);
%fprintf(1,'about to pause\n');
%pause;
