function [minerr,maxerr]=is_conj_sym(R)
%function [minerr,maxerr]=is_conj_sym(R)

ndim1=size(R,1);
ndim2=size(R,2);
testval=zeros(size(R));
for x1=1:ndim1 %index into Matlab matrix
  x2=conj_sym_index_1D(x1,ndim1);
  for y1=1:ndim2 %index into Matlab matrix
    y2=conj_sym_index_1D(y1,ndim2);
    if x1~=x2 | y1~=y2 %not going to count a possible imaginary part as error
      err=R(x1,y1)-conj(R(x2,y2));
    else
      if imag(R(x1,y1))~=0.0
        fprintf(1,'(%d,%d) is self conjugate symmetric but has value %g+i*%g\n', ...
                x1,y1,real(R(x1,y1)),imag(R(x1,y1)));
      end
      err=0.0;
    end


    denominator=abs(R(x1,y1))+abs(R(x2,y2));
    if denominator~=0.0
      testval(x1,y1)=abs(err)./denominator;
    elseif abs(err)==0.0
      testval(x1,y1)=0.0;
    else
      testval(x1,y1)=abs(err)./denominator; %really are dividing by 0
    end
  end
end
minerr=min(testval(:));
maxerr=max(testval(:));
