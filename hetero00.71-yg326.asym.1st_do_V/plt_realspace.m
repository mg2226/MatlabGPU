function [rhobar,rxx,x_rect]=plt_realspace(wantrhobar,wantrxx,onevobj,onetilde_b,m1l,m1h,d1,m2l,m2h,d2,m3l,m3h,d3,alpha,beta,gamma)
%function [rhobar,rxx,x_rect]=plt_realspace(wantrhobar,wantrxx,onevobj,onetilde_b,m1l,m1h,d1,m2l,m2h,d2,m3l,m3h,d3,alpha,beta,gamma)

%This function is almost able to also compute 3-D reciprocal space
%cubes because DoerschukJohnson IEEE T-IT 2000 Eq. 3 is so similar to
%Eq. 2.  The angular basis functions are the same, but I might change
%the name of the variable 'angular' to 'aAngular' to emphasize this
%point.  The radial basis functions must be changed: replace
%hlpr0_c_vec by hlpk0_c_vec and then multiply each basis function or
%each c_{l,n,p} value by (-1)^l.  This multiplication is not difficult
%because the basis functions and c_{l,n,p} are both real.  However,
%the answer will be complex so rhobar would be complex.  That is not a
%problem in matlab, but maybe I would do something else because of an
%interest in C.  For rxx there is more to think about because there
%are really multiple correlation functions depending on where you put
%\Re and \Im or ^\ast.  I have done similar calculations in the past,
%e.g., LeeDoerschukJohnson IEEE T-IP 2007 Eqs. 51-54.

%Set in case they are not computed.
rhobar=[];
rxx=[];
x_rect=[];

R2=onevobj.R2;
ll=onevobj.clnp.il;
mm=onevobj.clnp.in; %Spherical harmonics notation of 'm' not icosahedral harmonics notation of 'n'.
pp=onevobj.clnp.ip;
c=onevobj.clnp.c;
nu=onevobj.nu;
WhichAngularBasis=onevobj.BasisFunctionType;
tilde_b=onetilde_b.tilde_b;

SphericalBasis=0;
IcosahedralBasis=1;

%Initialize variables for hlpr0/hlpk0.
rtmax=10000.0;
lmax=max(ll);
pmax=max(pp);
if R2>rtmax
  error('plt_realspace: R2 %g > rtmax %g\n',R2,rtmax);
end
[root,aa]=init_hlp0_Hlp0_c(R2,lmax,pmax,rtmax);

%Setup the grid.
[m1,m2,m3]=ndgrid(m1l:m1h,m2l:m2h,m3l:m3h);
m123=[m1(:).*d1, m2(:).*d2, m3(:).*d3];
Rinv=euler2R([-gamma -beta -alpha]);
x_rect=m123*Rinv; %rotate the coordinates by Rinv
x_sph=xyztosph_vec(x_rect); %convert rectangular coordinates to spherical coordinates

[unique_lp,map_lp2unique,map_unique2lp]=unique([ll pp],'rows');
[unique_radius,map_radius2unique,map_unique2radius]=unique(x_sph(:,1));
htable=zeros(length(unique_radius),length(unique_lp));
for lp=1:length(unique_lp)
  htable(:,lp)=hlpr0_c_vec(unique_lp(lp,1),unique_lp(lp,2),unique_radius,root,aa,R2);
end

%rhobar: ZhengWangDoerschuk JOSA-A 2012 Eq. 16.
%rxx: ZhengWangDoerschuk JOSA-A 2012 Eq. 18 with $\vec x^\prime=\vec x$ and
%$(V_{\eta^\prime})_{\tau,\tau^\prime}=\nu_{\eta^\prime}(\tau)\delta_{\tau,\tau^\prime}$.

if wantrhobar
  rhobar=zeros(size(x_sph,1),1);
end
if wantrxx
  rxx=zeros(size(x_sph,1),1);
end

for ii=1:length(ll)

  l=ll(ii);
  m=mm(ii);

  %Not yet implemented: reuse h when l and p have not changed from the previous iteration.
  h=htable(:,map_unique2lp(ii));
  h_full=h(map_unique2radius);

  %Reuse angular if ii>1 && l==ll(ii-1)
  if ii==1 || l~=ll(ii-1)

    if WhichAngularBasis==IcosahedralBasis
      angular=set_Ttable_nord(l,m,x_sph(:,2),x_sph(:,3),tilde_b,onevobj.Ftable);
    elseif WhichAngularBasis==SphericalBasis
      angular=set_Ytable(l,x_sph(:,2),x_sph(:,3));
    else
      error('plt_realspace: WhichAngularBasis %d #1\n',WhichAngularBasis);
    end

  end 

  if WhichAngularBasis==IcosahedralBasis
    if wantrhobar
      rhobar=rhobar + c(ii) .* ( (h_full .* angular{1}(:))    );
    end
    if wantrxx
      rxx=rxx + nu(ii) .*      ( (h_full .* angular{1}(:)).^2 );
    end
  elseif WhichAngularBasis==SphericalBasis
    if wantrhobar
      rhobar=rhobar + c(ii) .* ( (h_full .* angular{1}(mm(ii)+1,:)')    );
    end
    if wantrxx
      rxx=rxx + nu(ii) .*      ( (h_full .* angular{1}(mm(ii)+1,:)').^2 );
    end
  else
    error('plt_realspace: WhichAngularBasis %d #2\n',WhichAngularBasis);
  end

end %close for ii=1:length(ll)

if ~isempty(rhobar)
  rhobar=reshape(rhobar,m1h-m1l+1,m2h-m2l+1,m3h-m3l+1);
end
if ~isempty(rxx)
  rxx=reshape(rxx,m1h-m1l+1,m2h-m2l+1,m3h-m3l+1);
end
