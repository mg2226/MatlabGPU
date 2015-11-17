function Alndet = lndet(A)
%function Alogdet = lndet(A)
%
% Yili Zheng
% 8/31/2007
%

% pivoting in LU factorization may result in negative diagonal entries.
% so we need to take the absolute value of the diagnal.
[AL, AU] = lu(A);
Alndet = sum(log(abs(diag(AU)))); 