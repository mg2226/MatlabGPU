function Ytable=set_Ytable(ll,theta,phi)
%function Ytable=set_Ytable(ll,theta,phi)

Nc=length(ll);
Ny=length(theta);
assert (length(phi)==Ny);
Ytable = cell(Nc,1);
costheta=cos(theta);
sqrt2 = sqrt(2);

for ii=1:Nc
  l = ll(ii);
  if ii>1 && l==ll(ii-1)
    Ytable{ii} = Ytable{ii-1}; 
  else
    plgndrtable = legendre(l,costheta);
    % Ylm is stored in the m+1 index because the first index in Matlab is 1
    % but m starts with 0.
    Ytable{ii} = zeros(2*l+1,Ny);
    % m=0 case
    Nl0 = sqrt((2*l+1)/(4*pi));
    Ytable{ii}(1,:) = Nl0*plgndrtable(1,:);
    for m=2:2:2*l
      Nlm = sqrt((2*l+1)*factorial(l-m/2)/(4*pi*factorial(l+m/2)));
      m_2_phi = m/2.*phi(:);
      Ytable{ii}(m,:) = sqrt2*Nlm*plgndrtable(m/2+1,:)'.*sin(m_2_phi); % m={1,3,...,2l-1}
      Ytable{ii}(m+1,:) = sqrt2*Nlm*plgndrtable(m/2+1,:)'.*cos(m_2_phi); % m={2,4,...,2l}
    end
  end
end