function [avg_pixel_value, roi] = getAvgPixValROI(image, window, roiposition, cord, sPath)
% EK 23

%% align the image data with contour+outline map using given coordinates

f = round(window/2)+1; % first stimulus evoked response frame 

load([sPath(1:end-1) filesep num2str(1) filesep 'contours.mat']) %
load([sPath(1:end-1) filesep num2str(1) filesep 'coordinates.mat'], 'snap')
contours = imread([sPath(1:end-1) filesep num2str(1) filesep 'overlaid_contourMap.png']);
contours = rgb2ind(contours, 256);

allData = imresize(image, [size(contours,1), size(contours,2)]); % resize the experiment image to the size of ret exp 
ref = allData(:,:,f);
fixedpts=cord.fixed;
movingpts=cord.moved;

if isempty(roiposition)

    [ref_moved,~] = imalign_s(ref, contours,movingpts, fixedpts);
    
    figure;
    imagesc(ref_moved); hold on; freezeColors
    contour(azi, azi_levels, 'b', 'LineWidth', 1); daspect([1 1 1]);  hold on;
    contour(nanmean(imout,3), 'k', 'LineWidth',0.2); axis ij; axis off
    axis image
    
   h = imfreehand;
   roi = createMask(h); % roi as a mask  
else 
    roi = roiposition; % use the given roi position
end 

%% Calculate the average pixel value within the ROI
avg_pixel_values = zeros(window, 1);

% Iterate over each frame
for i = 1:size(allData,3) 

[temp,~] = imalign_s(allData(:,:,i), contours, movingpts, fixedpts);
pixValues = temp(roi);
avg_pixel_value(i) = mean(pixValues(:));
end
end 
