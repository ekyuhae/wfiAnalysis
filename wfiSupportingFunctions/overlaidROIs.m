% Create a figure
figure;
a = [0.4 0.6 0.7]
% Plot the rectangular ROIs
% for i = 1:length(roi)
%     subplot(
%     rectangle('Position', [roi{i}(1), roi{i}(2), roi{i}(3), roi{i}(4)]);
%     
% end

%%
figure; %rois are cell in barmapping and structure in gratings

    % subplot(4,5,i);
     subplot(1,2,1);
    % Plot the data for the current subplot 
        imagesc(avgStimResponse.monoc{1,end}); hold on;
    % Get the current axis limits
    ax = gca;
    xlimits = ax.XLim;
    ylimits = ax.YLim;
    colormap viridis
    % Draw the ROI rectangle
    % rectangle('Position', roi.monoc, 'EdgeColor', 'r', 'LineWidth', 2);
    h = imshow(cat(3, roi.binoc, zeros(size(roi.binoc)), zeros(size(roi.binoc))));
set(h, 'AlphaData', 0.5); % Set transparency (adjust as needed)

    % Set the axis limits to match the original data
    xlim(xlimits);
    ylim(ylimits);
    % Set the title for each subplot
    title(sprintf('ROI %d', i));
    subplot(1,2,2)
imagesc(avgStimResponse.binoc{1,end});