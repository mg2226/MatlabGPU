function [g,a,b]=scale_rhobar(map1,map2,capsid)
% Calculate a linear scaling between two maps over the capsid region 
% Factor g minimizes sum(map1-g*map2)^2
% Factors a and b minimize sum(map1-a*map2-b)^2 

if size(map1)~=size(map2)
    error('Input maps do not have the same sizes');
end

map12=map1.*map2;
sum12=sum(sum(sum(map12(capsid))));
map2sq=map2.^2;
sum2sq=sum(sum(sum(map2sq(capsid))));
g=sum12/sum2sq;

count=sum(sum(sum(capsid)));
sum2=sum(sum(sum(map2(capsid))));
sum1=sum(sum(sum(map1(capsid))));
A=zeros(2,2);
A(1,1)=sum2;
A(1,2)=count;
A(2,1)=sum2sq;
A(2,2)=sum2;
c=[sum1
   sum12];
x=A\c;
a=x(1);
b=x(2);