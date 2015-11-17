function Ttable=set_Ttable_nord(ll, nn, theta, phi,tilde_b,Ftable)
%function Ttable=set_Ttable(ll, nn, theta, phi)
%
% Compute icosahedral harmonics basis functions
%
% See double ti(int l, int n, double phi, double theta) in ti.c
%
% Yili Zheng
% 06/25/2008
%
% b_fn = 'callti.out.read_by_C'; % file contains \tilde{b} values computed from mathematica 
% lmax = max(ll);
% tilde_b = rd_b(b_fn, lmax);

Nc=length(ll);
Ny=length(theta);
assert (length(phi)==Ny);
Ttable = cell(Nc,1);
costheta=cos(theta);

prev_l = -1;
prev_n = -1;
for ii=1:Nc
  l = ll(ii);
  n = nn(ii);
  if l==prev_l && n==prev_n
    Ttable{ii} = Ttable{ii-1}; 
  else
    %plgndrtable = legendre(l,costheta);
    % Tlm is stored in the m+1 index because the first index in Matlab is 1
    % but m starts with 0.
    Ttable{ii} = zeros(Ny,1);
    if (mod(l,2)==0) % l even
      for m=n:l/5
        mprime = m*5;
        Nlm = sqrt((2*l+1)*Ftable(l-mprime+1)/(4*pi*Ftable(l+mprime+1)));
        if m==0
          %Ttable{ii}(:) = Ttable{ii}(:) + Nlm.*tilde_b(l+1,n+1,m+1).*plgndrtable(mprime+1,:)'.*cos(mprime*phi); 
          Ttable{ii}(:) = Ttable{ii}(:) + Nlm.*tilde_b(l+1,n+1,m+1).*plgndr(l,mprime,costheta).*cos(mprime*phi);
%             Ttable{ii}(:,n) = Ttable{ii}(:,n) + Nlm.*tilde_b(l+1,n+1,m+1).*plgndr(l,mprime,costheta).*cos(mprime*phi);
        else
          %Ttable{ii}(:) = Ttable{ii}(:) + 2*Nlm.*tilde_b(l+1,n+1,m+1).*plgndrtable(mprime+1,:)'.*cos(mprime*phi); 
          Ttable{ii}(:) = Ttable{ii}(:) + 2*Nlm.*tilde_b(l+1,n+1,m+1).*plgndr(l,mprime,costheta).*cos(mprime*phi);
        end
      end
    else % l odd
      for m=n+1:l/5
        mprime = m*5;
        Nlm = sqrt((2*l+1)*Ftable(l-mprime+1)/(4*pi*Ftable(l+mprime+1)));
        %Ttable{ii}(:) = Ttable{ii}(:) + 2*Nlm*tilde_b(l+1,n+1,m+1)*plgndrtable(mprime+1,:)'.*sin(mprime*phi);
        Ttable{ii}(:) = Ttable{ii}(:) + 2*Nlm.*tilde_b(l+1,n+1,m+1).*plgndr(l,mprime,costheta).*sin(mprime*phi);
      end
    end
    prev_l = l;
    prev_n = n;
  end
end
