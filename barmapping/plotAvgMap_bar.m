function plotAvgMap_bar(contrasts, avgStimResponse, loc, coordinates, sPath)
% after pulling all sessions' data, get avg stim response datasets and plot
% them altogether
% EK 23


load([sPath(1:end-1) filesep num2str(1) filesep 'contours.mat'])
load([sPath(1:end-1) filesep num2str(1) filesep 'coordinates.mat'], 'snap')
contours = imread([sPath(1:end-1) filesep num2str(1) filesep 'overlaid_contourMap.png']);
contours = rgb2ind(contours, 256);

allData = imresize(avgStimResponse, [size(contours,1), size(contours,2)]);
for i = 1:size(allData,3)
    [avgStimResponse_resized(:,:,i),~] = imalign_s(allData(:,:,i), contours, coordinates.moved,coordinates.fixed);
end
%% 
m = figure;
% colorscale = [min(min(min(avgStimResponse))) max(max(max(avgStimResponse)))];
tiledlayout('flow','TileSpacing', 'compact', 'Padding', 'compact'); % show avg stim response activity for each white bar location
for i = 1: length(loc)
    nexttile
%     clims = [min(avgStimResponse.binoc(:)) max(avgStimResponse.binoc(:))];
    imagesc(avgStimResponse_resized(:,:,i)); 
    colormap viridis; colorbar; freezeColors
    title([num2str(loc(i)) ' degrees']);

    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);   %changed_AR
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off %changed_AR
    axis image
    
end   
sgtitle ([sprintf('%d%%', floor(contrasts)), ' contrast bars']); 

cd(sPath)

set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 25])
savefig(m, ['avg barmap activity_' num2str(contrasts) ' perc.fig']);
saveas(m, ['avg barmap activity_' num2str(contrasts) ' perc.png']);
