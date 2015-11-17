function L  =setL_nord(Rabc,ll,mm,vk,Htable,map_unique2lp,tilde_b)
%function L=setL(Rabc,ll,mm,vk,Htable,map_unique2lp)
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

%fprintf(1,'This is function setL\n');

Nc=length(ll);
Ny=size(vk,1);
L=zeros(Ny*2,Nc);


kk=[vk, zeros(Ny,1)]; % kk is the array of 3D vectors that includes the z corrdinate, which is 0
kknew=kk*Rabc; % apply rotation to all k vectors
% Rabc' is better than Rabc when used on synthetic images with test_l3p5
% /cacRuns/Test_asym_rule/
% fprintf(1,'Rabc is \n');
% disp(Rabc);
% kknew=kk*Rabc'; % apply rotation to all k vectors

sphkk=xyztosph_vec(kknew); % convert cartesian coordinates to spherical coordinates for all k vectors

% set up the table of spherical harmonics 
%unique_ll = unique(ll);
% Ytable=set_Ytable(ll,sphkk(:,2),sphkk(:,3));
Atable=set_Ttable_nord(ll,mm,sphkk(:,2),sphkk(:,3),tilde_b);
% Atable=set_Ytable(ll,sphkk(:,2),sphkk(:,3));

% set up L 
for ii=1:Nc
  switch (mod(ll(ii),4))
    case 0 % (-i)^0=1
%
        L(1:2:Ny*2-1,ii) = Atable{ii}(:).*Htable(:,map_unique2lp(ii)); %
%         real part for Ttable
%         L(1:2:Ny*2-1,ii) = Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))'; % real part for Ytable
     
    case 1 % (-i)^1=-i
%       
      L(2:2:Ny*2,ii) = -Atable{ii}(:).*Htable(:,map_unique2lp(ii)); % imaginary part
%         L(2:2:Ny*2,ii) = -Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))'; % imaginary part
      
    case 2 % (-i)^2=-1
%       
      L(1:2:Ny*2-1,ii) = -Atable{ii}(:).*Htable(:,map_unique2lp(ii)); % real part
%         L(1:2:Ny*2-1,ii) = -Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))'; % real part
      
    case 3 % (-i)^3=i
%       
      L(2:2:Ny*2,ii) = Atable{ii}(:).*Htable(:,map_unique2lp(ii)); % imaginary part
%         L(2:2:Ny*2,ii) = Atable{ii}(mm(ii)+1,:).*Htable(:,map_unique2lp(ii))'; % imaginary part
      
    otherwise
      error('Incorrect value of mod(ll(ii),4) %d!', mod(ll(ii),4));
  end
end

