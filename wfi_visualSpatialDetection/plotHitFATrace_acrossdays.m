close all; clear; clc
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK\project'))
addpath(genpath("Y:\haider\Code\behaviorAnalysis"));


% area = {[], 'LM_PM', 'RL_AM'};
area = {[]};
fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\';
period = 'days'; %'it can be 'animals'
RT = '';
% cd([fPath filesep 'stats' filesep 'behavior'])
% cd('Y:\haider\Data\analyzedData\EKK\WFI\M230831_1\behavior\LM_PM_compiled_days')
for i  = 1:length(area)
    hva = area{i};
    if isempty(hva)
        hit = load(['compiled ROI traces_' period '_hit' RT '.mat']);
        FA = load(['compiled ROI traces_' period '_FA' RT '.mat']);
    else
        hit = load(['compiled ROI traces_' period '_hit_' hva RT '.mat']);
        FA = load(['compiled ROI traces_' period '_FA_' hva RT '.mat']);
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
    % if ~isempty(RT)
        xlabel('relative time to lick onset (s)')
    % end
    setfig()

    subplot 122
    [tracem,stdm] = getAvgTraceAcrossSess(temp1(1,2), hit.idxm, 30, 0);
    if isempty (hva)
        title('Monocular V1');
    else
        title(hva(4:5))
    end
    % if ~isempty(RT)
        xlabel('relative time to lick onset (s)')
    % end

    % sgtitle (['Average Response over N = ' num2str(height(explist)) ' days'])
    yLim = [-2 7];
    % yLim = [min(min(peak_value))-max(max(std_peak)) max(max(peak_value))+max(max(std_peak))]; % Define the y-axis limit
    ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
    set(ax,'ylim',yLim);
    setfig()
    %%
    type = 'hit';
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