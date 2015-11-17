function unitq=rotation2q(R)
% Recover unit quaternion q=(w,x,y,z) based on given 3x3 rotation matrix

qxx=R(1,1);
qyy=R(2,2);
qzz=R(3,3);
qxy=R(1,2);
qxz=R(1,3);
qyx=R(2,1);
qyz=R(2,3);
qzx=R(3,1);
qzy=R(3,2);


t=qxx+qyy+qzz;
r=sqrt(1+t);
w=0.5*r;
x=abs(0.5*sqrt(1+qxx-qyy-qzz));
if qzy-qyz<0
    x=-x;
end
y=abs(0.5*sqrt(1-qxx+qyy-qzz));
if qxz-qzx<0
    y=-y;
end
z=abs(0.5*sqrt(1-qxx-qyy+qzz));
if qyx-qxy<0
    z=-z;
end

unitq=[w,x,y,z];