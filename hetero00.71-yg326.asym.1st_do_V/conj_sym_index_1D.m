function jj=conj_sym_index_1D(ii,N)
%function jj=conj_sym_index_1D(ii,N)
%ii is the index to a Matlab vector: ii\in\{1,\dots,N\}

if ii<1 | ii>N
  error(['ii ' num2str(ii) ' N ' num2str(N)]);
end

jj=-(ii-1)+N+1;
if ii==1
  jj=ii;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
