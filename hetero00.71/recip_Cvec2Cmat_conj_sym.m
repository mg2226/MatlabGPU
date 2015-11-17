function Y3=recip_Cvec2Cmat_conj_sym(Y2,Na,Nb,ixminimalset,iyminimalset)
%function Y3=recip_Cvec2Cmat_conj_sym(Y2,Na,Nb,ixminimalset,iyminimalset);

n=size(Y2,1);

Y3=NaN*ones(Na,Nb);
filled=zeros(Na,Nb);

for m=1:n

  %(x1,y1) was computed, (x2,y2) is the conjugate symmetric point
  x1=ixminimalset(m);
  x2=conj_sym_index_1D(x1,Na);
  y1=iyminimalset(m);
  y2=conj_sym_index_1D(y1,Nb);
%  fprintf(1,'%s: (x1,y1) (%d,%d) (x2,y2) (%d,%d)\n',mfilename,x1,y1,x2,y2);

  if filled(x1,y1)==1 || filled(x2,y2)==1
    error(['m ' num2str(m)]);
  end

  Y3(x1,y1)=Y2(m);
  if x1~=x2 || y1~=y2 %in case there is a small imaginary part
    Y3(x2,y2)=conj(Y2(m));
  elseif imag(Y3(x1,y1))~=0.0
    fprintf(1,'%s: m %d (%d,%d) is self conjugate symmetric but has value %g+i*%g\n', ...
            mfilename,m,x1,y1,real(Y3(x1,y1)),imag(Y3(x1,y1)));
    Y3(x1,y1)=real(Y3(x1,y1));
  end
  filled(x1,y1)=1;
  filled(x2,y2)=1;

end

minfilled=min(filled(:));
maxfilled=max(filled(:));
if minfilled~=1 || maxfilled~=1
  fprintf(1,'%s: min(filled(:)) %g max(filled(:)) %g\n', ...
          mfilename,minfilled,maxfilled);
end

[minerr,maxerr]=is_conj_sym(Y3);
if minerr~=0.0 || maxerr~=0.0
  fprintf(1,'%s: Y3 conjugate symmetry: minerr %g maxerr %g\n', ...
          mfilename,minerr,maxerr);
end
