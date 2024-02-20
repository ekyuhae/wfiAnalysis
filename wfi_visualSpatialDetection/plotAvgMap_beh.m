function plotAvgMap_beh (contrasts, avgStimResponse, coordinates, sPath, type, hva, reactionTriggered)
% after pulling all sessions' data, get avg stim response datasets and plot
% them altogether
% EK 23
% NA 6/1/2023 - added catch for hvas
% NA 6/2/2023 - changed plot titles to have hva, added hva as last input

%Used in: [1] behSpatioTempResp.m, [2] behSpatioTempResp_event.m
%%
sessDir = fileparts(pwd); 
% if ~strcmp(sessDir, [sessDir filesep num2str(1)]) && length(sessDir) > length([sessDir(1:end-2) '\' num2str(1)]) %if it's subfolder
%   sessDir = [sessDir '\..'];
% end

if ~isempty(hva) || reactionTriggered
    sessDir = [sessDir '\..'];
end
    load([sessDir filesep num2str(1) filesep 'contours.mat']) %
    load([sessDir filesep num2str(1) filesep 'coordinates.mat'], 'snap')
contours = imread([sessDir filesep num2str(1) filesep 'overlaid_contourMap.png']);
% catch
%     load([sPath filesep 'contours.mat'])
%     load([sPath filesep 'coordinates.mat'])
%     contours = imread([sPath filesep 'overlaid_contourMap.png']);
% % end
contours = rgb2ind(contours, 256);

%%
idx = find(contrasts == 0, 1, 'last');

bcont = contrasts(1:idx-1);
for i = 1: length(avgStimResponse.binoc)
    if ~isempty(avgStimResponse.binoc{i})
        binoc{i} = avgStimResponse.binoc{i};
    else
        binoc{i}= NaN;
    end
end

m = figure;
% colorscale = [min(min(min(binoc))) max(max(max(binoc)))];
% tiledlayout(1, length(contrasts),'TileSpacing', 'tight', 'Padding', 'compact'); % show avg stim response activity for each white bar location
tiledlayout(1, length(bcont),'TileSpacing', 'compact', 'Padding', 'compact'); % show avg stim response activity for each white bar location
for i = 1: length(bcont)
    binoc_resized{i} = imresize(binoc{i}, [size(contours,1), size(contours,2)]);
    [binoc_resized{i},~] = imalign_s(binoc_resized{i}, contours, coordinates.moved,coordinates.fixed);
    
    nexttile
    imagesc(binoc_resized{i}); colormap viridis; freezeColors; colorbar
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    
    title(sprintf('%d%%', bcont(i)*100));
end
if isempty(hva)
    sgtitle (['avg activity map - binoc ' type ' trials']);
else
    sgtitle(['avg activity map - ' hva(1:2) ' ' type ' trials'])
end

set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
if isempty(hva)
    savefig(m,  [sPath filesep 'avgActivityMap_binoc_' type 'Trials.fig']);
else 
    savefig(m, [sPath filesep 'avgActivityMap_' hva(1:2) '_' type 'Trials.fig'])
end
% saveas(m,  ['avgActivityMap_binoc_' type 'Trials.png']);


mcont = contrasts(idx:end);

for i = 1: length(avgStimResponse.monoc)
    if ~isempty(avgStimResponse.monoc{i})
        monoc{i} = avgStimResponse.monoc{i};
    else
        monoc{i}= NaN;
    end
end

f = figure;
% colorscale = [min(min(min(monoc))) max(max(max(monoc)))];
% tiledlayout(1, length(contrasts),'TileSpacing', 'tight', 'Padding', 'compact'); % show avg stim response activity for each white bar location
tiledlayout(1, length(mcont),'TileSpacing', 'compact', 'Padding', 'compact'); % show avg stim response activity for each white bar location
for i = 1: length(mcont)
    monoc_resized{i} = imresize(monoc{i}, [size(contours,1), size(contours,2)]);
    [monoc_resized{i},~] = imalign_s(monoc_resized{i}, contours, coordinates.moved,coordinates.fixed);
    
    nexttile
    imagesc(monoc_resized{i}); colormap viridis; freezeColors; colorbar
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    title(sprintf('%d%%', mcont(i)*100));
end
if isempty(hva)
    sgtitle (['avg activity map - monoc ' type ' trials']);
else
    sgtitle (['avg activity map - ' hva(4:5) ' ' type ' trials'])
end
colorbar
set(f, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
if isempty(hva)
    savefig(f, [sPath filesep 'avgActivityMap_monoc_' type 'Trials.fig']);
else
    savefig(f, [sPath filesep 'avgActivityMap_' hva(4:5) '_' type 'Trials.fig'])
end
% saveas(f,  ['avgActivityMap_monoc_' type 'Trials.png']);
