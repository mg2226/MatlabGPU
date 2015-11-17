function capsid=select_capsid(map1,map2,r,thres)
% Determine the region of capsid based on two rhobar maps. Select all the
% voxels constrained by inner radius r(1) and outerradius r(2) and have
% intensity above certain percentage of the maximum.
% Output is a binary matrix 

s1=size(map1);
s2=size(map2);
if s1~=s2
    error('Two maps do not have the same sizes');
end
center=round(s1/2);
max1=max(max(max(map1)));
max2=max(max(max(map2)));
for i=1:s1(1)
    for j=1:s1(2)
        for k=1:s1(3)
            capsid(i,j,k)=false;
            dis=sqrt((i-center(1))^2+(j-center(2))^2+(k-center(3))^2);
            if dis>=r(1) && dis<=r(2)
                if (map1(i,j,k)>=thres(1)*max1)&&(map2(i,j,k)>=thres(2)*max2)
                    capsid(i,j,k)=true;
                end
            end
        end
    end
end

% The threshold and radii for capsid may need to be manully adjusted based 
% on the cross section plots.
crossxy=zeros(s1(1:2));
crossxy(:,:)=capsid(:,:,center(3));
figure;
imshow(crossxy);

crossyz=zeros(s1(2:3));
crossyz(:,:)=capsid(center(1),:,:);
figure;
imshow(crossyz);

crossxz=zeros(s1(1),s1(3));
crossxz(:,:)=capsid(:,center(2),:);
figure;
imshow(crossxz);

