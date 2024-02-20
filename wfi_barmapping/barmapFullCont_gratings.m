function [binocAvgStimResponse, monocAvgStimResponse, binocFrames, monocFrames] = barmapFullCont_gratings(paramSorted, allData, bFrameTimes, sPath, window)
% bar mapping full contrast analysis. updated version of WFIbarmappingAnalysis2
% plot activity maps corresponding to different bar locations
% EK Feb23
%% 
binocTime = paramSorted{3}(3:end);
monocTime = paramSorted{4}(3:end);
FATime_binoc = paramSorted{1}(3:end);
FATime_monoc = paramSorted{2}(3:end);

%% WFI barmapping analysis

% binocAvgStimResponse = zeros(size(allData,1), size(allData,2), size(binocTime,2));
% monocAvgStimResponse = zeros(size(allData,1), size(allData,2), size(monocTime,2));
%   blueStim = find((bTimes - stimOn) > 0, 1);

gratingsOn = [];
binocFrames = [];
monocFrames = [];

gratingsOn = findStimOnFrame(binocTime, 50, bFrameTimes); 
for i = 1:length(gratingsOn)
    binocFrames = cat(2,binocFrames, gratingsOn(i)-window:gratingsOn(i)+window);
end
binocAvgStimResponse = nanmean(allData(:, :, binocFrames), 3);


gratingsOn = findStimOnFrame(monocTime, 50, bFrameTimes); 
for i = 1:length(gratingsOn)
    monocFrames = cat(2,monocFrames, gratingsOn(i)-window:gratingsOn(i)+window);
end
monocAvgStimResponse = nanmean(allData(:, :, monocFrames), 3);
%%
m = figure; % show avg stim response activity for each white bar location
s = subplot (1,2,1);
clims = [min(binocAvgStimResponse(:)) max(binocAvgStimResponse(:))];
imagesc(binocAvgStimResponse, clims);  colormap viridis;
title (['avg response to binoc gratings']); %axis square;
colorbar
s1= subplot (1,2,2);
clims = [min(monocAvgStimResponse(:)) max(monocAvgStimResponse(:))];
imagesc(monocAvgStimResponse, clims);  colormap viridis
title (['avg response to monoc gratings']); %axis square;
colorbar
    
set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 30 10]) 
savefig(m, [sPath filesep 'gratings response.fig']);
saveas(m, [sPath filesep 'gratings response.png']);

%% 
% [m, idx] = max(binocAvgStimResponse(:,:,length(wTime)), [], "all", "linear"); %find max value of last frame bar position
% [pix1, pix2, ~] = ind2sub(size(binocAvgStimResponse),idx); %find corresponding pixel
% wpixelval = squeeze(nanmean(reshape(binocAvgStimResponse(pix1, pix2, :), [], length(wTime)),1)); %find the pixel value for each bar position
% 
% [m, idx] = max(monocAvgStimResponse(:,:,length(bTime)), [], "all", "linear");
% [pix1, pix2, ~] = ind2sub(size(binocAvgStimResponse),idx);
% bpixelval = squeeze(nanmean(reshape(monocAvgStimResponse(pix1, pix2, :), [], length(bTime)),1));
% 
% f = figure; 
% plot (paramSorted(2,1:length(wTime)), wpixelval, '-o'); hold on
% plot (paramSorted(2,1:length(bTime)), bpixelval, '-o')
% ylabel('dF/F'); xlabel('bar location')
% legend ('white', 'black')
% title ('Stimulus triggered activity at different bar locations')
% savefig([sPath filesep 'pixelValue.fig']);
% saveas(f,[sPath filesep 'pixelValue.png']);
