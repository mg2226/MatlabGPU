function [vavg,unique_radius]=spherical_avg_v(vobj,tilde_b)
%Calculate spherical averages of V based on Eq.37 in Q.Wang et
%al./JSB(2013)

ll=vobj.clnp.il;
pp=vobj.clnp.ip;
nn=vobj.clnp.in;
lmax = max(ll);
pmax = max(pp);
%ll, pp (1060,1) vectors
%tilde_b is a chart padded with zeros, size(56,2,12)
%tilde_b = rd_b('callti.out.read_by_C',lmax); %modified by pd83 08-07-2015: pass tilde_b rather than read tilde_b
[ul,map_ll2ul,map_ul2ll]=unique(ll);
nl=length(ul);% nl is the number of unique l values, not the same to the Nl in Eq.37 in Q.Wang et al./JSB(2013)
bln=zeros(nl,2);
nu=vobj.nu;
if length(nu)~=length(ll)
    fprintf('wrong size read nu');
end
% sum|blnm|^2
for il=1:nl
    l=ul(il);
    n=[0,1];
    m=0;
    tempb=[abs(tilde_b(l+1,n(1)+1,m+1)).^2,abs(tilde_b(l+1,n(2)+1,m+1)).^2];
    for m=1:size(tilde_b,3)-1 
        % b(l,n,m)=b(l,n,-m)*(-1)^(l+m)
        tempb=tempb+[(abs(tilde_b(l+1,1,m+1)).^2)*2,(abs(tilde_b(l+1,2,m+1)).^2)*2];
    end
    bln(il,:)=tempb;
end
% Htable; from plt_realspace.m
R2=280;
rtmax=10000.0;
if R2>rtmax
  error('plt_realspace: R2 %g > rtmax %g\n',R2,rtmax);
end
[root,aa]=init_hlp0_Hlp0_c(R2,lmax,pmax,rtmax);

%Setup the grid.
m1l=-R2;
m1h=R2;
m2l=-R2;
m2h=R2;
m3l=-R2;
m3h=R2;
d1=1;
d2=1;
d3=1;
gamma=0;
beta=0;
alpha=0;
[m1,m2,m3]=meshgrid(m1l:m1h,m2l:m2h,m3l:m3h);
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

%
indln2list=[];
[unique_ln,map_ln2unique,map_unique2ln]=unique([ll nn],'rows');
for i=1:size(unique_ln,1)
    if unique_ln(i,2)==1
        indln2list=[indln2list unique_ln(i,1)];
    end
end
indln2=unique(indln2list);
vlpb=zeros(length(unique_lp),1);
up=unique(pp);
for i=1:length(unique_lp)
    ib=ceil(i/length(up));
    templ=ul(ib);
    %indln2=[30 36 40 42 45];
    if ~isempty(find(templ==indln2))
        vlpb(i)=nu(map_lp2unique(i))*bln(ib,2);
        vlpb(i)=vlpb(i)+nu(map_lp2unique(i)-length(up))*bln(ib,1);
    else
        vlpb(i)=nu(map_lp2unique(i))*bln(ib,1);
    end
end

%
vavg=zeros(length(unique_radius),1);
for i=1:length(unique_lp)
        vavg(:)=vavg(:)+vlpb(i)*(htable(:,i).^2);
end
vavg=vavg/4/pi;
