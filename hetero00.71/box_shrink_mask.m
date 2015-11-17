function mask=box_shrink_mask(dimensions,pixels2delete)
%function mask=box_shrink_mask(dimensions,pixels2delete)
%dimensions(1:2) -- non-negative integers, dimensions of the image
%pixels2delete(1:4) -- non-negative integers, number of rows/columns of pixels to delete
%the order is "left", "bottom", "right", "top" when you think of the image as a matrix.

if any(pixels2delete<0)
  error('box_shrink_mask: pixels2delete %d %d %d %d\n',pixels2delete(1),pixels2delete(2), ...
        pixels2delete(3),pixels2delete(4));
end

mask=ones(dimensions);
mask(:,1:pixels2delete(1))=0;
mask(end-pixels2delete(2)+1:end,:)=0;
mask(:,end-pixels2delete(3)+1:end)=0;
mask(1:pixels2delete(4),:)=0;
mask=logical(mask);
