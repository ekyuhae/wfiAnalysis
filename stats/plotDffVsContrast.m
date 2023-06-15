function f = plotDffVsContrast(bcontrast, mcontrast, peak_value, std, expID)

f = figure;
subplot(1,2,1)
errorbar(bcontrast*100, peak_value(:,1), std(:,1), 'ko', 'LineWidth', 1.5);

xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('binocular'); axis square; axis padded; box off

subplot(1,2,2)
errorbar(mcontrast*100, peak_value(:,2), std(:,2), 'ko', 'LineWidth', 1.5);

xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('monocular'); axis square; axis padded; box off
if ~isempty(expID)
sgtitle(['average response activity across sesssions - ' expID]);
end 

yLim = [min(min(peak_value))-max(max(std)) max(max(peak_value))+max(max(std))]; % Define the y-axis limit
ax = findobj(f,'type','axes'); % Find all axes handles in the figure
set(ax,'ylim',yLim);

