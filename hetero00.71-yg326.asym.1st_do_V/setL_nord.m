function L=setL_nord(Rabc,ll,mm,vk,Htable,map_unique2lp,tilde_b,Ftable,BasisFunctionType)
%function L=setL_nord(Rabc,ll,mm,vk,Htable,map_unique2lp,tilde_b,Ftable,BasisFunctionType)
%
% Rabc: the rotation matrix of the euler angle (alpha, beta, gamma)
% ll: the array of all l values
% pp: the array of all p values
% mm: the array of all m values
% vk: the array of all kappa pairs
% respectively
%
% Yili Zheng
% 2/16/2008
%
% l and m start from 0!

% Corrected version
% Use Ttable to replace Ytabble
% Corrected Line 44 to Line 58 in computing L from Htable and Ytable
% pass in tildeb instead of reading it every time
%
% Qiu Wang
% 07/01/2010

%modified by pd83 on 07-07-2015 by allowing choice of angular basis functions:
%icosahedral harmonics or real-valued spherical harmonics.
%Icosahedral: use plgndr to compute P_{l,m}(\cos\theta) for scalar l, scalar m, and vector \theta.
%Spherical: use matlab's legendre to compute P_{l,m}(\cos\theta) for scalar l, all m, and vector \theta.
%This means that Atable{ii} is a vector for icosahedral but is a matrix for real-valued spherical.
%
%Ftable is set in virusobj_set_2Dreciprocal in all situations.
%
%The coefficients tilde_b are read only by the operator 'EM_read_tilde_b'.  The array can be empty if they are not required.

%fprintf(1,'This is function setL\n');

Nc=length(ll);
Ny=size(vk,1);
L=zeros(Ny*2,Nc);

kk=[vk, zeros(Ny,1)]; % kk is the array of 3D vectors that includes the z corrdinate, which is 0
kknew=kk*Rabc; % apply rotation to all k vectors
% Correction 08/2014 : use Rabc instead of Rabc' to avoid 5-repetition pattern in pdfs over
% orientations
% Rabc' is better than Rabc when used on synthetic images with test_l3p5
% /cacRuns/Test_asym_rule/
% fprintf(1,'Rabc is \n');
% disp(Rabc);
% kknew=kk*Rabc'; % apply rotation to all k vectors

%sphkk=xyztosph_vec(kknew); % convert cartesian coordinates to spherical coordinates for all k vectors
[sphkk(:,3),sphkk(:,2),sphkk(:,1)]=cart2sph(kknew(:,1),kknew(:,2),kknew(:,3));
sphkk(:,2)=pi/2-sphkk(:,2);

%set up the table of angular basis functions:
switch BasisFunctionType
  case 0 %real-valued spherical harmonics
    %I could first do
    %[unique_ll,unique2ll,ll2unique]=unique(ll);
    %which implies
    %ll(ii)=unique_ll(ll2unique(ii));
    %and
    %unique_ll(ii)=ll(unique2ll(ii));
    %Passing unique_ll will compute a smaller Atable.  But, so long as
    %the c_{l,n,p} are listed such that coefficients with the same
    %value of l are sequential, the absence of 'unique' makes the data
    %structure larger but doesn't cost a lot of time.  If I include
    %'unique', then in the four cases controlled by (-i)^l, I must
    %look up the correct Atable{??} to use through Atable{ll2unique(ii)} .
    Atable=set_Ytable(ll,sphkk(:,2),sphkk(:,3));
  case 1 %icosahedral harmonics
    Atable=set_Ttable_nord(ll,mm,sphkk(:,2),sphkk(:,3),tilde_b,Ftable);
  otherwise
    fprintf(1,'setL_nord: BasisFunctionType %d\n',BasisFunctionType);
    error('setL_nord: Unknown BasisFunctionType');
end

%set up L 
for ii=1:Nc
  switch (mod(ll(ii),4))

    case 0 % (-i)^0=1, therefore real part with a "+" sign
      switch BasisFunctionType
        case 0 %real-valued spherical harmonics
          L(1:2:Ny*2-1,ii) = Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))';
        case 1 %icosahedral harmonics
          L(1:2:Ny*2-1,ii) = Atable{ii}(:).*Htable(:,map_unique2lp(ii));
        otherwise
          fprintf(1,'setL_nord: BasisFunctionType %d\n',BasisFunctionType);
          error('setL_nord: Unknown BasisFunctionType');
      end
     
    case 1 % (-i)^1=-i, therefore imaginary part with a "-" sign
      switch BasisFunctionType
        case 0 %real-valued spherical harmonics
          L(2:2:Ny*2,ii) = -Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))';
        case 1 %icosahedral harmonics
          L(2:2:Ny*2,ii) = -Atable{ii}(:).*Htable(:,map_unique2lp(ii));
        otherwise
          fprintf(1,'setL_nord: BasisFunctionType %d\n',BasisFunctionType);
          error('setL_nord: Unknown BasisFunctionType');
      end
      
    case 2 % (-i)^2=-1, therefore real part with a "-" sign
      switch BasisFunctionType
        case 0 %real-valued spherical harmonics
          L(1:2:Ny*2-1,ii) = -Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))';
        case 1 %icosahedral harmonics
          L(1:2:Ny*2-1,ii) = -Atable{ii}(:).*Htable(:,map_unique2lp(ii));
        otherwise
          fprintf(1,'setL_nord: BasisFunctionType %d\n',BasisFunctionType);
          error('setL_nord: Unknown BasisFunctionType');
      end
      
    case 3 % (-i)^3=i, therefore imaginary part with a "+" sign
      switch BasisFunctionType
        case 0 %real-valued spherical harmonics
          L(2:2:Ny*2,ii) = Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))';
        case 1 %icosahedral harmonics
          L(2:2:Ny*2,ii) = Atable{ii}(:).*Htable(:,map_unique2lp(ii));
        otherwise
          fprintf(1,'setL_nord: BasisFunctionType %d\n',BasisFunctionType);
          error('setL_nord: Unknown BasisFunctionType');
      end

    otherwise
      error('Incorrect value of mod(ll(ii),4) %d!', mod(ll(ii),4));
  end
end
