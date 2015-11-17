function [tildeb, space_check, ll, Ipl]=rd_b_general(b_fn, lmax)
% Extend the "rd_b.m" function for the use of all 5 Icosahedral group classes (Ip=1,2,3,4,5)
%Input:
% b_fn: the file name of tildeb data file. The file should have the format
% below:
%     l1    //First l value
%     Ip1   //1st Ip value which has nonzero basis function coefficients under
%           //the above l1
%     {}    //1st nonzero basis coefficients vector under l1, Ip1
%     ...
%     {}    //last nonzero basis coefficients vector under l1, Ip1 
%     l1    //First l value
%     Ip2   //2nd Ip value which has nonzero basis function coefficients under
%           //l1
%     {}    //1st nonzero basis coefficients vector under l1, Ip2
%     ...
%     {}    //last nonzero basis coefficients vector under l1, Ip2 
%     l2    //Second l value
%     Ip1'  //1st Ip value which has nonzero basis function coefficients under
%           //l2
%     {}    //1st nonzero basis coefficients vector under l2, Ip1'
%     ...
%     {}    //last nonzero basis coefficients vector under l2, Ip1'
%
% lmax: max l

%Output: 
% tildeb: a cell array of size (lmax+1) X Ip, where each cell entry
% includes a matrix of basis coefficient matrix. Row of the 
% matrix is the vector of basis coefficient for m=-l,...,l.
%
% space_check: a (lmax+1) X 1 vector to test if the total number of basis
% function under certain l span the whole space: (0, yes; 1 no.)
%
% ll: a vector of list of all ls 
% Ipl: a vector of list of all Ips
%
% Test:
% b_fn='RealBasisFunctionCoeff_RealUniIrrepPhyPiover2_noStructure_1012.txt'; lmax=45;
% Nan Xu
% 10/16/2014


nIp=5;
tildeb=cell(lmax+1, nIp); %lmax+1-->l, nIp-->p, nmax-->given l, each Ip has n basis, 2l
space_check=zeros(lmax+1,1);
[fid,fopenmsg]=fopen(b_fn,'r');
if fid == -1
  error(['rd_b_general: fid ' num2str(fid) ' fopenmsg ' fopenmsg ' opening ' b_fn]);
end

TotalLength=2*sum(0:lmax)+lmax+1;
ll=zeros(TotalLength,1); Ipl=zeros(TotalLength,1); 
tline = fgetl(fid); rd_l1=str2num(tline); 
tline = fgetl(fid); rd_Ip=str2num(tline);  %read the l and Ip
if isempty(rd_Ip) || isempty(rd_Ip)
    error('rd_b_general: file opening reading error\n');
end
% if (rd_l1~=0 || rd_Ip~=1)
%     error('rd_b_general: l 0 rd_l %d Ip 1 rd_Ip %d\n',rd_l1,rd_Ip);
% end
I=1i; ll_ct=0;
for l=0:lmax     
    n=0;
    rd_l2=l;
    if (rd_l1~=l || rd_Ip>5)
        error('rd_b_general: l %d rd_l %d Ip 1~5 rd_Ip %d\n',l, rd_l1,rd_Ip);
    end
    while rd_l1==rd_l2 && ~feof(fid)
      tline = fgetl(fid); %read the 1st basis function coefficient vector      
      tildeb_curm=[];
      while strcmp(tline(1),'{')==1
        ll_ct=ll_ct+1; ll(ll_ct)=l; Ipl(ll_ct)=rd_Ip;
        tildeb_curv=cell2mat(eval(tline));
        tildeb_curm=[tildeb_curm; tildeb_curv];
        tline = fgetl(fid); 
      end
      tildeb{l+1,rd_Ip}=tildeb_curm;
      n=n+size(tildeb_curm,1);
      if  ~feof(fid)
          rd_l2=str2num(tline); %tline is the new l
          tline = fgetl(fid);      
          rd_Ip=str2num(tline); %tline is the new Ip
      end
    end
    rd_l1=rd_l2;
    
    if n~=2*l+1
        space_check(l+1)=1;
    end
end
fclose(fid);
