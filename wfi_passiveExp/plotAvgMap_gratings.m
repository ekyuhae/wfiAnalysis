function plotAvgMap_gratings (contrasts, avgStimResponse, coordinates, sPath, hva)
% after pulling all sessions' data, get avg stim response datasets and plot
% them altogether
% EK 23, NA 7/11/23 added hva
% used in: barmapMultiCont_gratings.m
%% 

load([sPath(1:end-1) filesep num2str(1) filesep 'contours.mat'])
load([sPath(1:end-1) filesep num2str(1) filesep 'coordinates.mat'], 'snap')
contours = imread([sPath(1:end-1) filesep num2str(1) filesep 'overlaid_contourMap.png']);
contours = rgb2ind(contours, 256);

%%
for i = 1: length(avgStimResponse.binoc)
    if ~isempty(avgStimResponse.binoc{i})
        binoc(:,:,i) = avgStimResponse.binoc{i};
    else
        binoc(:,:,i)= NaN;
    end
end

for i = 1: length(avgStimResponse.monoc)
    if ~isempty(avgStimResponse.monoc{i})
        monoc(:,:,i) = avgStimResponse.monoc{i};
    else
        monoc(:,:,i)= NaN;
    end
end

m = figure;
% colorscale = [min(min(min(binoc))) max(max(max(binoc)))];
% tiledlayout(1, length(contrasts),'TileSpacing', 'tight', 'Padding', 'compact'); % show avg stim response activity for each white bar location
tiledlayout(1, length(contrasts),'TileSpacing', 'compact', 'Padding', 'compact'); % show avg stim response activity for each white bar location
for i = 1: length(contrasts)
    binoc_resized(:,:,i) = imresize(binoc(:,:,i), [size(contours,1), size(contours,2)]);
    [binoc_resized(:,:,i),~] = imalign_s(binoc_resized(:,:,i), contours, coordinates.moved,coordinates.fixed);
    
    nexttile
    imagesc(binoc_resized(:,:,i)); colormap viridis; freezeColors; colorbar    
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);   
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off 
    axis image
    
    title(sprintf('%d%%', contrasts(i)*100));    
end 
sgtitle ('avg response to binoc gratings'); 

set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(m, ['avg response to binoc gratings' hva '.fig']);
saveas(m,  ['avg response to binoc gratings' hva '.png']);

f = figure; 
% colorscale = [min(min(min(monoc))) max(max(max(monoc)))];
% tiledlayout(1, length(contrasts),'TileSpacing', 'tight', 'Padding', 'compact'); % show avg stim response activity for each white bar location
tiledlayout(1, length(contrasts),'TileSpacing', 'compact', 'Padding', 'compact'); % show avg stim response activity for each white bar location
for i = 1: length(contrasts)
    monoc_resized(:,:,i) = imresize(monoc(:,:,i), [size(contours,1), size(contours,2)]);
    [monoc_resized(:,:,i),~] = imalign_s(monoc_resized(:,:,i), contours, coordinates.moved,coordinates.fixed);
    
    nexttile
    imagesc(monoc_resized(:,:,i)); colormap viridis; freezeColors; colorbar
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    title(sprintf('%d%%', contrasts(i)*100));
end
sgtitle ('avg response to monoc gratings'); 
colorbar
set(f, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(f, ['avg response to monoc gratings' hva '.fig']);
saveas(f, ['avg response to monoc gratings' hva '.png']);