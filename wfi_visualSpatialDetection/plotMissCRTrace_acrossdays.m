close all; clear; clc
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK\project'))


area = {'RL_AM'};
fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\';
period = 'days'; %'it can be 'animals'

% Animal = 'M230220_1';
RT = '';
% cd([fPath filesep 'stats' filesep 'behavior'])
for i  = 1:length(area)
    hva = area{i};
    if isempty(hva)
        hit = load(['compiled ROI traces_days_miss' RT '.mat']);
        FA = load(['compiled ROI traces_days_CR' RT '.mat']);
    else
        hit = load(['compiled ROI traces_days_miss_' hva RT '.mat']);
        FA = load(['compiled ROI traces_days_CR_' hva RT '.mat']);
    end
    %%
    for i =1:2
        temp1{1,i} = zeros(size(hit.avg_exp{1,1}));
        if isequal(size(hit.avg_exp{1,1}), size(FA.avg_exp{1,1}))
            temp1{1,i}(1,:,:) = FA.avg_exp{1,i}(1,:,:);
            temp1{1,i}(2:end,:,:) = hit.avg_exp{1,i}(2:end,:,:);
        else
            temp1{1,i}=zeros(size(hit.avg_exp{1,1}));
            x = cat(3, FA.avg_exp{1,i}(1,:,:), repmat(NaN,1,91,1));
            temp1{1,i}(1,:,:) = x(1,:,:);
            temp1{1,i}(2:end,:,:) = hit.avg_exp{1,i}(2:end,:,:);
        end
        % temp1{1,i}(2:end,:,:) = hit.temp{1,i}(2:end,:,:);
    end

    f1 = figure;

    subplot 121
    [traceb,stdb] = getAvgTraceAcrossSess(temp1(1,1), hit.idxb, 30, 0);
    if isempty (hva)
        title('Binocular V1');
    else
        title(hva(1:2))
    end
    if ~isempty(RT)
        xlabel('relative time to lick onset (s)')
    end
    setfig()

    subplot 122
    [tracem,stdm] = getAvgTraceAcrossSess(temp1(1,2), hit.idxm, 30, 0);
    if isempty (hva)
        title('Monocular V1');
    else
        title(hva(4:5))
    end
    if ~isempty(RT)
        xlabel('relative time to lick onset (s)')
    end

    % sgtitle (['Average Response over N = ' num2str(height(explist)) ' days'])
    yLim = [-2 7];
    % yLim = [min(min(peak_value))-max(max(std_peak)) max(max(peak_value))+max(max(std_peak))]; % Define the y-axis limit
    ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
    set(ax,'ylim',yLim);
    setfig()
    
    % for i =1:2
    %     temp1{1,i} = zeros(size(hit.temp{1,1}));
    %     temp1{1,i}(1,:,:) = FA.temp{1,i}(1,:,:);
    %     temp1{1,i}(2:end,:,:) = hit.temp{1,i}(2:end,:,:);
    % end
    % 
    % f1 = figure;
    % 
    % subplot 121
    % [traceb,stdb] = getAvgTraceAcrossSess(temp1(1,1), hit.data(1).idxb, 30, 0);
    % if isempty (hva)
    %     title('Binocular V1');
    % else
    %     title(hva(1:2))
    % end
    % if ~isempty(RT)
    %     xlabel('relative time to lick onset (s)')
    % end
    % 
    % subplot 122
    % [tracem,stdm] = getAvgTraceAcrossSess(temp1(1,2), hit.data(1).idxm, 30, 0);
    % if isempty (hva)
    %     title('Monocular V1');
    % else
    %     title(hva(4:5))
    % end
    % if ~isempty(RT)
    %     xlabel('relative time to lick onset (s)')
    % end
    % 
    % % sgtitle (['Average Response over N = ' num2str(height(explist)) ' days'])
    % yLim = [-1.5 4];
    % % yLim = [min(min(peak_value))-max(max(std_peak)) max(max(peak_value))+max(max(std_peak))]; % Define the y-axis limit
    % ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
    % set(ax,'ylim',yLim);
    %%
    type = 'miss';
    if isempty(hva)
        % cd([directory filesep 'V1_compiled_days'])
        savefig(f1, ['meandffTrace_animals_' type RT '.fig'])
        saveas(f1, ['meandffTrace_animals_' type RT '.png'])
    else
        % cd([directory filesep hva '_compiled_days'])
        savefig(f1, ['meandffTrace_animals_' type '_' hva RT '.fig'])
        saveas(f1, ['meandffTrace_animals_' type '_' hva RT '.png'])
    end

    clear temp1
end