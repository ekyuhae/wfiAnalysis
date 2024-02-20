clear; close all; clc

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
% Animal = 'M230717_1';
% period = 'days';
% date = '09-Oct-2023';
% attnblk = 'B';
% % flag = '';
% % flag = '_roiUnmatched';
% flag = '_blockUnmatched';
% if ~strcmp(period, 'days')
%     fPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior' filesep date];
% else
%     fPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior'];
% end
fPath = pwd;
hva = 'RL_AM';
fig1 = openfig(['meandff_over_contrasts_days_hit_' hva '.fig']);
fig2 = openfig(['meandff_over_contrasts_days_FA_' hva '.fig']);

axes1 = findobj(fig1, 'Type', 'axes');
data1 = get(axes1, 'Children'); % Get plot objects
for i = 1:2
    x_data1{i} = get(data1{i}, 'XData'); % Get x data
    y_data1{i} = get(data1{i}, 'YData'); % Get y data
    y_delta1{i} = get(data1{i}, 'YPositiveDelta');
end

% Figure 2 FA
axes2 = findobj(fig2, 'Type', 'axes');
data2 = get(axes2, 'Children'); % Get plot objects
for i = 1:2
    x_data2{i} = get(data2{i}, 'XData'); % Get x data
    y_data2{i} = get(data2{i}, 'YData'); % Get y data
    y_delta2{i} = get(data2{i}, 'YPositiveDelta');
end

% axes3 = findobj(fig3, 'Type', 'axes');
% data3 = get(axes3, 'Children'); % Get plot objects
% x_data3 = get(data3, 'XData'); % Get x data
% y_data3 = get(data3, 'YData'); % Get y data

%%
close all
% Plot the extracted data points from both figures together
figure;
for i = 1:2
    subplot(1,2,i)
    errorbar(x_data1{i}(2:end), y_data1{i}(2:end), y_delta1{i}(2:end),'k', 'lineWidth', 2, 'CapSize',0); % Plot data from figure 1
    hold on;
    errorbar(x_data1{i}(1), y_data1{i}(1), y_delta2{i}(1), 'lineWidth', 2,'CapSize',0); % Plot data from figure 1
    % xlim([0 10])
    ylim([-2 8])

    xlabel('Contrast(%)');
    ylabel('dF/F(%)');
    if i == 1
        if isempty (hva)
            title('Binocular V1');
        else
            title(hva(1:2))
        end
    else
        if isempty (hva)
            title('Monocular V1');
        else
            title(hva(4:5))
        end
    end
setfig()
    axis square; legend boxoff; box off; axis padded
    
    legend({'hit', 'FA'}, "Location",'southeast')
end

setfig()
cd(fPath)
savefig(['meandff_contrasts_hitFA.fig']);



%%
% close all
% % Plot the extracted data points from both figures together
% figure;
% 
% plot(x_data1{1}, y_data1{1}, 'color', [0.4 0.3 0.2], 'lineWidth', 2); % Plot data from figure 1
% hold on;
% plot(x_data2{1}, y_data2{1}, 'color', [0.7 0.1 0.3],  'lineWidth', 2); % Plot data from figure 1
% plot(x_data3{1}, y_data3{1}, 'color', [0.3 0.7 0.2], 'lineWidth', 2); % Plot data from figure 1
% 
% %baseline
% plot(x_data1{2}, y_data1{2}, 'color', [0.4 0.3 0.2], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
% hold on;
% plot(x_data2{2}, y_data2{2}, 'color', [0.7 0.1 0.3], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
% plot(x_data3{2}, y_data3{2}, 'color', [0.3 0.7 0.2], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
% 
% % poststim
% plot(x_data1{3}, y_data1{3},  'color', [0.4 0.3 0.2],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
% hold on;
% plot(x_data2{3}, y_data2{3},  'color', [0.7 0.1 0.3],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
% plot(x_data3{3}, y_data3{3},  'color', [0.3 0.7 0.2],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
% 
% xlim([0 15])
% ylim([-0.03 0.08])
% 
% xlabel('Hit trial number');
% ylabel('dF/F');
% if ~strcmp(flag, '_blockSwitched')
%     title('Post and Prestim dF/F in Monoc Attentional Sessions');
% else
%     title('Post and Prestim dF/F in the other block');
% end
% legend({'V1', 'LM', 'RL'});
% axis square; legend boxoff; box off
% setfig();
% 
% cd(fPath)
% savefig(['attentionAcrossAreas_Monoc' flag '.fig'])
