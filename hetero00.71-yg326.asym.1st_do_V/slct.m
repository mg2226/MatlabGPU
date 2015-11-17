function rt=slct(l,n)
%function rt=slct(l,n)
% slct.m modified frofunction rt=slct(l,n)m slct.c
% /* 2/21/96
%    Selection rule of icosahedral harmonics - if l-30*n is a non negative 
%    integer other than 1,2,3,4,5,7,8,9,11,13,14,17,19,23,29 then it is a 
%    valid order. 
%    Example: slct(0,0)=1, slct(6,0)=1, slct(30,1)=1, slct(3,0)=0.
% */

if (l - n * 30 < 0 || n < 0)
  rt = 0;
else
  switch (l - n * 30)
    case 1
      rt =  0;
    case 2
      rt =  0;
    case 3
      rt =  0;
    case 4
      rt =  0;
    case 5
      rt =  0;
    case 7
      rt =  0;
    case 8
      rt =  0;
    case 9
      rt =  0;
    case 11
      rt =  0;
    case 13
      rt =  0;
    case 14
      rt =  0;
    case 17
      rt =  0;
    case 19
      rt =  0;
    case 23
      rt =  0;
    case 29
      rt =  0;
      
    otherwise
      rt = 1;
      
  end
end