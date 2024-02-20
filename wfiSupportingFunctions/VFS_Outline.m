function VFS_Outline(VFS, azimuth, elevation, Animal, sPath)
%%creating outline of WFI VFS map 
% EK 22

addpath("\\ad.gatech.edu\bme\labs\haider\Code\Shared Code\isiAnalysisTonyNew\supporting_functions")
gradmag = abs(VFS); 
threshSeg = 2.0*std(VFS(:)); % thresholding VFS 
imseg = (sign(gradmag-threshSeg/2) + 1)/2; 
id = find(imseg);
imdum = imseg.*VFS; imdum(id) = imdum(id)+1.1;
patchSign = getPatchSign(imseg,VFS);
id = find(patchSign ~= 0);
patchSign(id) = sign(patchSign(id) - 1);
SE = strel('disk',2,0);
imseg = imopen(imseg,SE);
patchSign = getPatchSign(imseg,VFS);
imout = ploteccmap(patchSign,[1.1 2.1]);
    
h = figure;
% VFS with area overlay
subplot(1,3,1); imagesc(VFS, [-1 1]); colormap jet; colorbar; daspect([1 1 1]); hold on; 
contour(nanmean(imout,3), 'k');%, 'LineWidth',2);
title([Animal ' - VFS Contour Map']); axis image
% Azi with area overlay
subplot(1,3,2); imagesc(azimuth, [-30 110]);colormap jet;colorbar; daspect([1 1 1]); hold on;
contour(nanmean(imout,3),'k'); hold off; 
title([Animal ' - Azimuth+Contour']); axis image
% El with area overlay
subplot(1,3,3); imagesc(elevation, [-30 50]);colormap jet;colorbar; daspect([1 1 1]); hold on;
contour(nanmean(imout,3),'k'); hold off; 
title([Animal ' - Elevation+Contour']); axis image

savefig(h,[sPath filesep Animal '_contourMap.fig']);
saveas(h,[sPath filesep Animal '_contourMap.jpg']);
end 