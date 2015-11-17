function tildeb=rd_b(b_fn, lmax)
%function tildeb=rd_b(b_fn, lmax)
%modified from rd_b.c
%Input:
% b_fn: the file name of tildeb data file: 
% lmax: max l
%
% Example usage:
% tildeb=rd_b('/home/hours/yz247/research/EM4.4.ti.bigger_l.mathematica/callti.out.read_by_C',45)
%
% Yili Zheng
% Feb. 6, 2008
%
%
% Corrected by adding fclose()
%
% Qiu Wang
% 07/01/2010

nmax=floor(lmax/30);
mmax=floor(lmax/5);

tildeb=zeros(lmax+1,nmax+1,mmax+1);

[fid,fopenmsg]=fopen(b_fn,'r');
if fid == -1
  error(['rd_b: fid ' num2str(fid) ' fopenmsg ' fopenmsg ' opening ' b_fn]);
end

for l=0:lmax
    for n=0:l/30
      [tmp, tmp_count] = fscanf(fid,'%d%d',2);
      if ( tmp_count ~= 2)
        error('rd_b: l %d n %d reading l,n\n',l,n);
      end
      rd_l = tmp(1);
      rd_n = tmp(2);
      
      if (l~=rd_l || n~=rd_n)
        error('rd_b: l %d rd_l %d n %d rd_n %d\n',l,rd_l,n,rd_n);
      end

      if (slct(l,n)==1)
        if (mod(l,2)==0) % /* l even */
          for m=n:l/5
            [tmp, tmp_count] = fscanf(fid,'%lg',1);
            if ( tmp_count ~= 1)
              error('rd_b: l even, l %d n %d m %d\n',l,n,m);
            end
            tildeb(l+1,n+1,m+1)=tmp;
          end
        else % /* l odd */
          for m=n+1:l/5
            [tmp, tmp_count] = fscanf(fid,'%lg',1);
            if ( tmp_count ~= 1)
              error('rd_b: l odd, l %d n %d m %d\n',l,n,m);
            end
            % /*b[l][n+1][m]=tmp;*/
            tildeb(l+1,n+1,m+1)=tmp;
          end
        end
      end
    end
end
fclose(fid);
