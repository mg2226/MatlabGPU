function [LT_SigmaInv_y LT_SigmaInv_L]=LT_SigmaInv_fast(q,L,nu,y)
%function y_Lc_Sigma=y_Lc_Sigma_fast(q,L,nu,y_Lc)

% Qiu Wang
% 06/15/2010

Nc=length(nu);
%Nl=size(L,1);

if q==0
    Sigma=L*diag(nu)*L';
    LT_SigmaInv_y=L'/Sigma*y;
    LT_SigmaInv_L=L'/Sigma*L;
    return
elseif q<0
    fprintf(1,'Error: noise_var not positive \n');
    return
end

LTy=L'*y;
LTL=L'*L;
LT_SigmaInv_y=1/q*LTy;
LT_SigmaInv_L=1/q*LTL;
if nnz(nu)==Nc
    midmx=1/q*(L'*L);
    for i=1:Nc
        midmx(i,i)=midmx(i,i)+1/nu(i);
    end
    LT_SigmaInv_y=LT_SigmaInv_y-q^(-2)* (LTL/midmx)*LTy;
    LT_SigmaInv_L=LT_SigmaInv_L-q^(-2)* LTL* (midmx\LTL);
end