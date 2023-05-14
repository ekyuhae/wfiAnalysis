function VFS_Outline_Contour(VFS, azi, Animal, sPath)
%%creating outline of WFI VFS map + contour map of elevation and azimuth
%%retinotopy map
% EK 22

addpath("Y:\haider\Code\Shared Code\isiAnalysisTonyNew\supporting_functions")
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
    

%% 
% Load your azimuth and elevation gradient maps into the variables 'azimuth' and 'elevation'

% Define the contour levels you want to plot (in degrees)
azi_levels = -60:10:140;
% elev_levels = -40:10:50;

figure;

% Compute the contour lines for each level
contour(azi, azi_levels, 'LineWidth', 1);  hold on;
contour(nanmean(imout,3), 'k', 'LineWidth',0.2); axis ij; axis off; axis image
f = getframe;
colorbar; daspect([1 1 1]); 
f.cdata = imresize(f.cdata, [size(azi,1) size(azi,2)]);

save([sPath filesep 'contours.mat'], 'imout', 'azi', 'azi_levels');
imwrite(f.cdata,[sPath filesep  'overlaid_contourMap.png']);
savefig([sPath filesep 'overlaid_contourMap.fig']);
saveas(gcf, [sPath filesep 'overlaid_contourMap.jpg']); % image w contour level bar
end 