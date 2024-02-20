clear; close all; clc

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
Animal = 'M230831_1';
period = 'sessions';
date = '25-Oct-2023';
attnblk = 'M';
% flag = '';
% flag = '_ROIswitched';
flag = '_blockSwitched';
if ~strcmp(period, 'days')
    fPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior' filesep date];
else
    fPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior'];
end

if strcmp(attnblk, 'B')
    if ~strcmp(flag, '_blockSwitched') 
        cd([fPath filesep 'attention_over_' period])
        % fig1 = openfig('peak_vs_BinocHitTrialNumber.fig');
        fig1 = openfig(['peak_vs_BinocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_RL_AM'])
        % fig2 = openfig('peak_vs_BinocHitTrialNumber.fig');
        fig2 = openfig(['peak_vs_BinocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_LM_PM'])
        % fig3 = openfig('peak_vs_BinocHitTrialNumber.fig');
        fig3 = openfig(['peak_vs_BinocHitTrialNumber' flag '.fig']);

    else
        cd([fPath filesep 'attention_over_' period])
        fig1 = openfig(['peak_vs_MonocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_RL_AM'])
        fig2 = openfig(['peak_vs_MonocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_LM_PM'])
        fig3 = openfig(['peak_vs_MonocHitTrialNumber' flag '.fig']);

    end

    axes1 = findobj(fig1, 'Type', 'axes');
    data1 = get(axes1, 'Children'); % Get plot objects
    x_data1 = get(data1, 'XData'); % Get x data
    y_data1 = get(data1, 'YData'); % Get y data

    % Figure 2
    axes2 = findobj(fig2, 'Type', 'axes');
    data2 = get(axes2, 'Children'); % Get plot objects
    x_data2 = get(data2, 'XData'); % Get x data
    y_data2 = get(data2, 'YData'); % Get y data

    axes3 = findobj(fig3, 'Type', 'axes');
    data3 = get(axes3, 'Children'); % Get plot objects
    x_data3 = get(data3, 'XData'); % Get x data
    y_data3 = get(data3, 'YData'); % Get y data

    %%
    close all
    % Plot the extracted data points from both figures together
    figure;

    plot(x_data1{1}, y_data1{1}, 'color', [0.4 0.3 0.2], 'lineWidth', 2); % Plot data from figure 1
    hold on;
    plot(x_data2{1}, y_data2{1}, 'color', [0.7 0.1 0.3],  'lineWidth', 2); % Plot data from figure 1
    plot(x_data3{1}, y_data3{1}, 'color', [0.3 0.7 0.2], 'lineWidth', 2); % Plot data from figure 1

    %baseline
    plot(x_data1{2}, y_data1{2}, 'color', [0.4 0.3 0.2], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    hold on;
    plot(x_data2{2}, y_data2{2}, 'color', [0.7 0.1 0.3], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    plot(x_data3{2}, y_data3{2}, 'color', [0.3 0.7 0.2], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1

    % poststim
    plot(x_data1{3}, y_data1{3},  'color', [0.4 0.3 0.2],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    hold on;
    plot(x_data2{3}, y_data2{3},  'color', [0.7 0.1 0.3],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    plot(x_data3{3}, y_data3{3},  'color', [0.3 0.7 0.2],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1

    xlim([0 10])
    ylim([-0.03 0.08])

    xlabel('Hit trial number');
    ylabel('dF/F');
    if  ~strcmp(flag, '_blockSwitched')
        title('Post and Prestim dF/F in Binoc Attentional Sessions');
    else
        title('Post and Prestim dF/F in the other block');
    end
    legend({'V1', 'RL', 'LM'});
    axis square; legend boxoff; box off
    setfig();

    cd(fPath)
    savefig(['attentionAcrossAreas_Binoc' flag '.fig']);


else

    if ~strcmp(flag, '_blockSwitched')%|| ~strcmp(flag, '_blockUnmatched')
        cd([fPath filesep 'attention_over_' period])
        fig1 = openfig(['peak_vs_MonocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_RL_AM'])
        fig2 = openfig(['peak_vs_MonocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_LM_PM'])
        fig3 = openfig(['peak_vs_MonocHitTrialNumber' flag '.fig']);
    else
        cd([fPath filesep 'attention_over_' period])
        fig1 = openfig(['peak_vs_BinocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_RL_AM'])
        fig2 = openfig(['peak_vs_BinocHitTrialNumber' flag '.fig']);

        cd([fPath filesep 'attention_over_' period '_LM_PM'])
        fig3 = openfig(['peak_vs_BinocHitTrialNumber' flag '.fig']);
    end

    axes1 = findobj(fig1, 'Type', 'axes');
    data1 = get(axes1, 'Children'); % Get plot objects
    x_data1 = get(data1, 'XData'); % Get x data
    y_data1 = get(data1, 'YData'); % Get y data

    % Figure 2
    axes2 = findobj(fig2, 'Type', 'axes');
    data2 = get(axes2, 'Children'); % Get plot objects
    x_data2 = get(data2, 'XData'); % Get x data
    y_data2 = get(data2, 'YData'); % Get y data

    axes3 = findobj(fig3, 'Type', 'axes');
    data3 = get(axes3, 'Children'); % Get plot objects
    x_data3 = get(data3, 'XData'); % Get x data
    y_data3 = get(data3, 'YData'); % Get y data

    %%
    close all
    % Plot the extracted data points from both figures together
    figure;

    plot(x_data1{1}, y_data1{1}, 'color', [0.4 0.3 0.2], 'lineWidth', 2); % Plot data from figure 1
    hold on;
    plot(x_data2{1}, y_data2{1}, 'color', [0.7 0.1 0.3],  'lineWidth', 2); % Plot data from figure 1
    plot(x_data3{1}, y_data3{1}, 'color', [0.3 0.7 0.2], 'lineWidth', 2); % Plot data from figure 1

    %baseline
    plot(x_data1{2}, y_data1{2}, 'color', [0.4 0.3 0.2], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    hold on;
    plot(x_data2{2}, y_data2{2}, 'color', [0.7 0.1 0.3], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    plot(x_data3{2}, y_data3{2}, 'color', [0.3 0.7 0.2], 'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1

    % poststim
    plot(x_data1{3}, y_data1{3},  'color', [0.4 0.3 0.2],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    hold on;
    plot(x_data2{3}, y_data2{3},  'color', [0.7 0.1 0.3],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1
    plot(x_data3{3}, y_data3{3},  'color', [0.3 0.7 0.2],'lineStyle', ':', 'lineWidth', 1.5); % Plot data from figure 1

    xlim([0 15])
    ylim([-0.03 0.08])

    xlabel('Hit trial number');
    ylabel('dF/F');
    if   ~strcmp(flag, '_blockSwitched')
             title('Post and Prestim dF/F in Monoc Attentional Sessions');
    else 
            title('Post and Prestim dF/F in the other block');
    end 
    legend({'V1', 'AM', 'PM'});
    axis square; legend boxoff; box off
    setfig();

    cd(fPath)
    savefig(['attentionAcrossAreas_Monoc' flag '.fig'])
end