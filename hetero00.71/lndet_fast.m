function lndetSigma=lndet_fast(q,L,nu,lndetV)
%function lndetSigma=lndet_fast(q,L,nu,lndetV)

Nc=length(nu);
Ny=size(L,1);

if q==0.0
  lndetSigma=lndet(L*diag(nu)*L');
  assert(isfinite(lndetSigma),'Assertion lndet_fast #1');
  return;
elseif q<0.0
  error('lndet_fast: q %g is <0.0\n',q);
end

if all(nu>0.0)
  %standard case
  M=(L'*L)./q;
  for i=1:Nc
    M(i,i)=M(i,i)+1.0/nu(i);
  end
  lndetSigma=Ny*log(q)+lndetV+lndet(M);
else
  %assume all(nu==0.0) is true and therefore drop the L V L^T term and return \ln\det(Q)
  lndetSigma=Ny*log(q);
end
assert(isfinite(lndetSigma),'Assertion lndet_fast #2');

% Sylvester's determinant theorem
% \det(A + UW^T) = \det(A) \det(I + W^T A^{-1} U)
% Definition of \Sigma (ZhengWangDoerschuk JOSA-A 2012 Eq. 7):
% \Sigma = L V L^T + Q
% Definition of V:
% V=\diag(\nu)
% Definition of Q:
% Q=q I_{N_y}
% Note that V and Q are both symmetric and both full rank.
% The goal is to compute \det\Sigma and eventually \ln\det\Sigma.
% Match the formula for \Sigma with the theorem:
% A=Q
% U=L
% W^T=VL^T
% Then
% \det\Sigma = \det(Q) \det(I + VL^T Q^{-1} L)
% = \det(Q) \det(V[ V^{-1} + L^T Q^{-1} L])
% = \det(Q) \det(V) \det(V^{-1} + L^T Q^{-1} L)
% = q^{N_y} \det(V) \det(V^{-1} + L^T (1/q) I_{N_y} L)
% = q^{N_y} \det(V) \det(V^{-1} + (1/q) L^T L)
% = q^{N_y} \det(V) \det(\diag(1/\nu) + (1/q) L^T L)
% I could write
% \det(V) = \det(\diag(\nu)) = \prod_{i=1}^{N_c} \nu
% and eventually
% \ln\det(V) = \ln\det(\diag(\nu)) = \sum_{i=1}^{N_c} \ln(\nu)
% but \nu is constant over the abscissias so I will compute \ln\det(V)
% outside of this function and pass the result to this function.
% Therefore,
% \ln\det\Sigma
% = \ln(q^{N_y}) + \ln\det(V) + \ln\det(\diag(1/\nu) + (1/q) L^T L)
% = N_y \ln(q) + \ln\det(V) + \ln\det(\diag(1/\nu) + (1/q) L^T L)
% Define M by
% M = \diag(1/\nu) + (1/q) L^T L
% Then
% \ln\det\Sigma
% = N_y \ln(q) + \ln\det(V) + \ln\det(M)
