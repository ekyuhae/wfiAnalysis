function plotTempTrace_bar(avgTempTrace, contrasts, fr, location, multi)
% plot mean temporal trace across the trials within particular roi
% applied to barmapping exp 
% EK 23 

%% Set the blue color values
k = linspace(0, 0.8, length(contrasts));
fr = fr/2;
t = 0:1/fr:(length(avgTempTrace{1,1})-1)/fr; % convert frame into timescale 
t = t - t(fr*2+1); 
% t = t - t(floor(fr*0.5)+1); 

f = figure; tiledlayout('flow','TileSpacing', 'compact', 'Padding', 'compact');
if multi
    % show avg stim response activity for each white bar location
    for locc = 1: length(location.loc) % plot traces with different contrast for every location
        nexttile       
        for i = 1:size(avgTempTrace,1) % repeat for contrasts
            avgTempTrace{end,locc} = avgTempTrace{end,1}; % 0 contrast
            legend_names{i} = sprintf('%d%%', contrasts(i));
            color = [k(i), k(i), k(i)];
            plot(t, avgTempTrace{i,locc}, 'Color', color);
            xlim([t(15) t(end)])
            hold on
            axis square;
        end
        title([num2str(location.loc(locc)) ' degrees']);
        xline(t(fr*2+1), ':k'); hold off       
    end
    legend(legend_names); legend('boxoff');
else
    for locc = 1: length(location.loc) % plot traces with different contrast for every location
        nexttile        
        plot(t, avgTempTrace{1,locc}, 'Color', [k(1) k(1) k(1)]);
        xlim([t(15) t(end)])
        hold on
        plot(t, avgTempTrace{2,1}, 'Color', [k(2) k(2) k(2)]);
        axis square;
        title([num2str(location.loc(locc)) ' degrees']);
        xline(t(fr*2+1), ':k'); hold off
    end
    legend('100%', '0%'); legend('boxoff'); 
end


a = reshape(avgTempTrace, [],1);
nonEmptyC = a(~cellfun('isempty', a)); % Select only non-empty cells
if ~isempty(nonEmptyC)
    mins = min(cellfun(@min, nonEmptyC));
    maxs = max(cellfun(@max, nonEmptyC));
else
    mins = NaN;
    maxs = NaN;
end
    

yLim = [mins maxs];
ax = findobj(f,'type','axes'); % Find all axes handles in the figure
set(ax,'ylim',yLim);

sgtitle ('mean ROI temporal trace'); 
set(f, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 25])
savefig(f, 'mean ROI trace.fig')
saveas(f, 'mean ROI trace.png')
