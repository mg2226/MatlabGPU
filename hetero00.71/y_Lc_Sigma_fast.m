function y_Lc_Sigma=y_Lc_Sigma_fast(q,L,nu,y_Lc)
%function y_Lc_Sigma=y_Lc_Sigma_fast(q,L,nu,y_Lc)

% Qiu Wang
% 06/15/2010

Nc=length(nu);
%Nl=size(L,1);

if q==0
    Sigma=L*diag(nu)*L';
    y_Lc_Sigma=sum( y_Lc .* (Sigma\y_Lc) );
    return
elseif q<0
    fprintf(1,'Error: noise_var not positive \n');
    return
end

y_Lc_Sigma= 1/q*sum(y_Lc.^2);
if nnz(nu)== Nc
    midmx=1/q*L'*L;
    for i=1:Nc
        midmx(i,i)=midmx(i,i)+1/nu(i);
    end
    LTy=L'*y_Lc;

    y_Lc_Sigma=y_Lc_Sigma-q^(-2)*sum(LTy.*(midmx\LTy));
end