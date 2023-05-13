function plotTempTrace_beh(avgTempTrace, contrasts, fr)
% plot mean temporal trace across the trials within particular roi
% applied to grating stimulus exp 
% EK 23 

%% Set the blue color values
idx = find(contrasts == 0, 1, 'last');
bcont = contrasts(1:idx-1);
mcont = contrasts(idx:end);

k = linspace(0.8, 0, length(bcont));
fr = fr/2;
t = 0:1/fr:(length(avgTempTrace.binoc)-1)/fr; % convert frame into timescale 
t = t - t(fr*2+1); 

b= figure;
legend_names = cell(1, length(bcont));
for i = 1:length(bcont)
    legend_names{i} = sprintf('%d%%', bcont(i)*100);
    color = [k(i), k(i), k(i)];    
%     y_adj = imadjust(avgTempTrace.binoc(i,:), [], [], i/length(contrasts));
%     plot(t, y_adj, 'Color', color);
     plot(t, avgTempTrace.binoc(i,:), 'Color', color);
    hold on;

end
xline(t(fr*2+1), ':k') 
title( 'binoc ROI trace');
xlabel('time (s)'); ylabel('dF/F_0')
legend(legend_names); legend('boxoff'); 
yLim = [min(min(avgTempTrace.binoc)) max(max(avgTempTrace.binoc))];
ax = findobj(b,'type','axes'); % Find all axes handles in the figure
set(ax,'ylim',yLim);
hold off


k = linspace(0.8, 0, length(mcont));
m = figure;
legend_names = cell(1, length(mcont));
for i = 1:length(mcont)
    legend_names{i} = sprintf('%d%%', mcont(i)*100);
    color = [k(i), k(i), k(i)];
    plot(t, avgTempTrace.monoc(i,:), 'Color', color);
    hold on; 

end
xline(t(fr*2+1), ':k')
title( 'monoc ROI trace');
xlabel('time (s)'); ylabel('dF/F_0')
legend(legend_names);  legend('boxoff'); 

yLim = [min(min(avgTempTrace.monoc)) max(max(avgTempTrace.monoc))];
ax = findobj(m,'type','axes'); % Find all axes handles in the figure
set(ax,'ylim',yLim);

savefig(b, 'binoc ROI trace in multi-contrast.fig')
saveas(b, 'binoc ROI trace in multi-contrast.png')
savefig(m, 'monoc ROI trace in multi-contrast.fig')
saveas(m, 'monoc ROI trace in multi-contrast.png')

end