function zeroone=zero_one_cube_from_angularbasisfunc(deltatheta,deltaphi,Nxyz,p,j,l,n,b_fn,lmax_in_b)
%function zeroone=zero_one_cube_from_angularbasisfunc(deltatheta,deltaphi,Nxyz,p,j,l,n,b_fn,lmax_in_b)
%test with
%b_fn='RealBasisFunctionCoeff_RealUniIrrepPhyPiover2_noStructure_1012.txt';
%lmax_in_b=45;

%read \cal D
[tildeb,space_check,ll,Ipl]=rd_b_general(b_fn,lmax_in_b);
%get the particular basis function that was requested--use j as the row index and ignore n
calDlp=tildeb{l+1,p};
calDlpj=calDlp(j,:); %vector of \cal D coefficients for -l<=m<=+l
if size(calDlpj,1)~=1 || size(calDlpj,2)~=2*l+1
  error(['size(calDlpj) ' num2str(size(calDlj))]);
end

%will compute I_{p;j,l,n}(\theta,\phi) on the following angle values:
thetavalues=[0:deltatheta:pi];
phivalues=[0:deltaphi:2*pi];

%need Y_{l,m}(\theta,\phi) for
%the one fixed value of l,
%all values of m (-l<=m<=+l),
%and all values of (\theta,\phi)

%Matlab legendre function:
%only has 0<=m<=+l not -l<=m<=+l
%Matlab writes P_n^m while PCD writes P_{l,m}:
%Matlab m = PCD m
%Matlab n = PCD l
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
%   Y fixed l,phi varying m,theta
    Y(m+1,:)=sqrt(tmp*factorialvalues(l-m+1)/factorialvalues(l+m+1)) ...
             * P(m+1,:)*expimphi(m+1);
  end
% I fixed l,phi varying theta
% add the terms for m^\prime=0,...,l for all values of theta
  Ipjln(nphi,:)=calDlpj(l+1:end)*Y;
% add the terms for m^\prime=-l,...,-1 (note reversed order of calDlpj elements)
  signfromm=1-2*mod(1:l,2);
  Ipjln(nphi,:)=Ipjln(nphi,:)+(signfromm.*calDlpj(l:-1:1))*conj(Y(2:end,:));
end
minIpjln=min(Ipjln(:));
maxIpjln=max(Ipjln(:));
%make Ipjln range from 0 to 1
Ipjln=(Ipjln-minIpjln)/(maxIpjln-minIpjln);
%make Ipjln range from 0.5 to 1
Ipjln=0.5*(1+Ipjln);

figure;
mesh(thetavalues/(2*pi),phivalues/(2*pi),Ipjln) %thetavalues=columns of Ipjln, phivalues=rwos of Ipjln
title('Ipjln after scaling');
xlabel('\theta/(2\pi)'); %thetavalues=columns of Ipjln
ylabel('\phi/(2\pi)'); %phivalues=rows of Ipjln

fprintf(1,'about to pause\n');
pause;

oversize=1.05;
deltaxyz=2*oversize/Nxyz;
samples=[-oversize:deltaxyz:oversize];
[x,y,z]=meshgrid(samples,samples,samples);
%[azimuth,elevation,r]=cart2sph(X,Y,Z);
[phi,theta,radius]=cart2sph(x,y,z);
theta=pi/2-theta;

Ipjln_interp2=interp2(thetavalues,phivalues,Ipjln,theta(:),phi(:));
Ipjln2=reshape(Ipjln_interp2,length(samples),length(samples),length(samples));

zeroone=(radius<Ipjln2);

nz=floor(length(samples)/2);

figure;
mesh(samples,samples,Ipjln2(:,:,nz));
title(['Ipjln2: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);

fprintf(1,'about to pause\n');
pause;

figure;
imagesc(zeroone(:,:,nz));
axis equal
title(['zeroone: p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);

fprintf(1,'about to pause\n');
pause;

figure;
mesh(samples,samples,radius(:,:,nz));
title(['p ' num2str(p) ' j ' num2str(j) ' l ' num2str(l) ' n ' num2str(n) ' x-y slice at nz ' num2str(nz) ' samples(nz) ' num2str(samples(nz))]);

fprintf(1,'about to pause\n');
pause;
