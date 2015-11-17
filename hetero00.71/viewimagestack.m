for ii=1:length(imagestack)
  imagesc(imagestack{ii});
  colormap(gray);
  axis equal;
  axis off;
  pause;
end
