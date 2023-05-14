function [mean_dFFtrace, roiposition] = getAvgTemporalResponse (paddedResponseFrames, img, allData, roiposition, cord, sPath) 
% take padded frames that correspond to stimulus evoked response and get
% temporal average across the number of stimulus repeats
% EK 2023 

% window = floor(img.sRate*0.25)+1; % if padded as -2s~2s window including
% stim onset frame; 
window = img.sRate*2+1; % 8s window, frame rate = 30
% window = 15;
temp = reshape(paddedResponseFrames', window, [])'; %matrix to get n (number of repeats) traces 

idx = find(temp > length(allData));
if ~isempty(idx)
    temp(end,:) = [];
end 

% for i = 1:size(temp,2) % avg across number of repeats
%     avgTemp(:,:,i) = mean(allData(:,:,temp(:,i)),3);
% end

for i = 1:size(temp,2) % avg across number of repeats
     idx = size(allData,3) > temp(:,i);
    avgTemp(:,:,i) = nanmean(allData(:,:,temp(idx,i)),3);
end

if isempty (roiposition) % get mean trace
    [mean_dFFtrace, roiposition] = getAvgPixValROI(avgTemp, window, [], cord, sPath); % for first contrast, set roi 
else
    [mean_dFFtrace,roiposition] = getAvgPixValROI (avgTemp, window, roiposition, cord, sPath);
end
end