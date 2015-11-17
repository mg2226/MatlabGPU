function ifnotexistmkdir(dn)
%function ifnoexistmkdir(dn)
%
% Create the directory dn if it doesn't exist
%
% Yili Zheng
% 05/13/2008
%

if ~exist(dn,'dir')
  [SUCCESS,MESSAGE,MESSAGEID] = mkdir(dn);
  if ~SUCCESS
    error('Failed to create the dir %s, error message %s, id %d', ...
      dn, MESSAGE, MESSAGEID);
  end
end
