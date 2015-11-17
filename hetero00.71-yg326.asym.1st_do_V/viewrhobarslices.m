for ii=1:size(rhobar,3)
  imagesc(rhobar(:,:,ii));
  colormap(gray);
  axis equal;
  axis off;
  fprintf(1,'rhobar(:,:,%d)\n',ii);
  pause;
end

for ii=1:size(rxx,3)
  imagesc(rxx(:,:,ii));
  colormap(gray);
  axis equal;
  axis off;
  fprintf(1,'rxx(:,:,%d)\n',ii);
  pause;
end
