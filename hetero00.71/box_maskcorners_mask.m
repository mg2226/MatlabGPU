function mask=box_maskcorners_mask(dimensions,delta,radius)
%function mask=box_maskcorners_mask(dimensions,delta,radius)
%dimensions(1:2) -- non-negative integers, dimensions of the image
%delta(1:2) -- positive reals, sampling intervals
%radius -- positive real, outside of the radius is "the corners"

[X1,X2]=ndgrid(delta(1)*([0:dimensions(1)-1]-(dimensions(1)-1)/2), ...
               delta(2)*([0:dimensions(2)-1]-(dimensions(2)-1)/2));
R=sqrt(X1.^2 + X2.^2);
mask=(radius<=R);
