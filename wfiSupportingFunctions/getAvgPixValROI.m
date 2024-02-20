function [avg_pixel_value, roi] = getAvgPixValROI(image, fwindow, roiposition, cord, sPath, flag)
% EK 23
% NA 6/1/2023 - added catch for hvas

%% align the image data with contour+outline map using given coordinates

f = round(fwindow/2)+1; % first stimulus evoked response frame 
cd(sPath) % either RT subfolders, HVA subfolders, or just V1 
sessDir = fileparts(pwd); % end at session number if it was subfolder
foldparts = strsplit(sessDir, filesep);

if (length(foldparts{end}) <= 2 && ~strcmp(foldparts{end},'1')) || (strcmp(foldparts{end},'1')) %if sessDir does not end at the exp date and end at session number = it was a subfolder
    sessDir = [sessDir '\..']; % shoud end at experiment date
end
% elseif  length(sessDir) > length([sessDir filesep num2str(1)]) % if its hva


load([sessDir filesep num2str(1) filesep 'contours.mat']) %
load([sessDir filesep num2str(1) filesep 'coordinates.mat'], 'snap')
contours = imread([sessDir filesep num2str(1) filesep 'overlaid_contourMap.png']);
% catch
%     load([sPath filesep 'contours.mat'])
%     load([sPath filesep 'coordinates.mat'])
%     contours = imread([sPath filesep 'overlaid_contourMap.png']);
% end
contours = rgb2ind(contours, 256);

allData = imresize(image, [size(contours,1), size(contours,2)]); % resize the experiment image to the size of ret exp 
ref = allData(:,:,f);
fixedpts=cord.fixed;
movingpts=cord.moved;

if isempty(roiposition)
    disp(['choose ' flag ' ROI'])
    [ref_moved,~] = imalign_s(ref, contours,movingpts, fixedpts);
    
    figure;
    imagesc(ref_moved); hold on; freezeColors
    contour(azi, azi_levels, 'b', 'LineWidth', 1); daspect([1 1 1]);  hold on;
    contour(nanmean(imout,3), 'k', 'LineWidth',0.2); axis ij; axis off
    axis image
    
   h = drawfreehand; % changed from imfreehand 5/12/23
   roi = createMask(h); % roi as a mask  
else 
    roi = roiposition; % use the given roi position
end 

%% Calculate the average pixel value within the ROI
avg_pixel_value = zeros(fwindow, 1);

% Iterate over each frame
for i = 1:size(allData,3) 

[temp,~] = imalign_s(allData(:,:,i), contours, movingpts, fixedpts);
pixValues = temp(roi);
avg_pixel_value(i) = mean(pixValues(:));
end
end 
