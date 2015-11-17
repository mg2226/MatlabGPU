function plt_realspace_crossections(cube,function4title)
%function plt_realspace_crossections(cube,function4title)

figure;
ii=floor(size(cube,1)/2);
imagesc(squeeze(cube(ii,:,:)));
colorbar;
title([function4title '(' num2str(ii) ',:,:)']);

figure;
ii=floor(size(cube,2)/2);
imagesc(squeeze(cube(:,ii,:)));
colorbar;
title([function4title '(:,' num2str(ii) ',:)']);

figure;
ii=floor(size(cube,3)/2);
imagesc(squeeze(cube(:,:,ii)));
colorbar;
title([function4title '(:,:,' num2str(ii) ')']);
