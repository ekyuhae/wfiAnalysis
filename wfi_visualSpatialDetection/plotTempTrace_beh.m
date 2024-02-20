function plotTempTrace_beh(avgTempTrace, contrasts, fr, type, hva)
% plot mean temporal trace across the trials within particular roi
% applied to grating stimulus exp 
% EK 23 
% NA 6/2/2023 - added hva input/change to subplot titling

%Used in: [1] behSpatioTempResp.m, [2] behSpatioTempResp_event.m

%%
idx = find(contrasts == 0, 1, 'last');
bcont = contrasts(1:idx-1);
mcont = contrasts(idx:end);

k = linspace(0.8, 0, length(bcont));
fr = fr/2;
t = 0:1/fr:(length(avgTempTrace.binoc)-1)/fr; % convert frame into timescale 
t = t - t(fr*2+1); 

b= figure;
subplot (1,2,1)
legend_names = cell(1, length(bcont));
for i = 1:length(bcont)
    legend_names{i} = sprintf('%d%%', bcont(i)*100);
    color = [k(i), k(i), k(i)];    
%     y_adj = imadjust(avgTempTrace.binoc(i,:), [], [], i/length(contrasts));
%     plot(t, y_adj, 'Color', color);
     plot(t, avgTempTrace.binoc(i,:), 'Color', color);
    hold on;

end
xline(t(fr*2+1), ':k') ;
if isempty(hva)
    title(['binoc ROI trace - ' type ' trials']);
else
    title([hva(1:2) ' ROI trace - ' type ' trials'])
end
xlabel('time (s)'); ylabel('dF/F_0')
legend(legend_names); legend('boxoff'); axis square
% if avgTempTrace.binoc ~= 0
%     yLim = [min(min(avgTempTrace.binoc)) max(max(avgTempTrace.binoc))];
%     ax = findobj(b,'type','axes'); % Find all axes handles in the figure
%     set(ax,'ylim',yLim);
% end
hold off

subplot (1,2,2)
k = linspace(0.8, 0, length(mcont));
legend_names = cell(1, length(mcont));
for i = 1:length(mcont)
    legend_names{i} = sprintf('%d%%', mcont(i)*100);
    color = [k(i), k(i), k(i)];
    plot(t, avgTempTrace.monoc(i,:), 'Color', color);
    hold on; 

end
xline(t(fr*2+1), ':k');
if isempty(hva)
    title(['monoc ROI trace - ' type ' trials']);
else
    title([hva(4:5) ' ROI trace - ' type ' trials'])
end
xlabel('time (s)'); ylabel('dF/F_0')
legend(legend_names);  legend('boxoff'); axis square
% if avgTempTrace.monoc ~= 0
%     yLim = [min(min(avgTempTrace.monoc)) max(max(avgTempTrace.monoc))];
%     ax = findobj(b,'type','axes'); % Find all axes handles in the figure
%     set(ax,'ylim',yLim);
% end

%%
if ~isempty(type)
    savefig(b, ['ROI traces_' type ' trials.fig'])
%     saveas(b, ['ROI traces_' type ' trials.png'])
%     savefig(m, ['monoc ROI trace_' type ' trials.fig'])
%     saveas(m, ['monoc ROI trace_' type ' trials.png'])
else
    savefig(b, 'ROI traces_whole trials.fig')
%     saveas(b, 'ROI traces_whole trials.png')
%     savefig(m, 'monoc ROI trace_whole trials.fig')
%     saveas(m, 'monoc ROI trace_whole trials.png')

end