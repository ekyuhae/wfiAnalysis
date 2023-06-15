%% load individual session datasets that were analyzed from lv2 and compile them 
% modified version of stat_new by AR
% EK 23 

close all; clear all; clc;
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
% addpath(genpath('\\ad.gatech.edu\bme\labs\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\';
Animal = {'M230220_1', 'M230209_1'};
% expID = {'03-Feb-2023', '08-Feb-2023_2','18-Feb-2023', '09-Mar-2023'};
expID = {'09-May-2023', '11-May-2023', '12-May-2023_1'};
gratings = 2; % false if barmapping
if gratings == 2 
    type = 'hit';
    hva = 'LM_PM';
end 
% cont = [100 75 50 25 5 0]; % barmapping contrast (shortcut)
% cont = [0,0.02,0.05,0.1,0.33,0.5,0.65,0.8,1]*100; % grating cont
%%
for m = 1: length(Animal)
    dir = [fPath Animal{m} filesep 'behavior'];

    if isempty(hva)
        cd([dir filesep 'V1_compiled_days'])
        data(m) = load(['compiled ROI traces_days_' type '.mat']);
    else
        cd([dir filesep hva '_compiled_days'])
        data(m) = load(['compiled ROI traces_days_' type '_' hva '.mat']);
    end
end
if isempty(hva)
    load([dir filesep data(m).expID{end} filesep 'compiled analyzed data.mat']); % load last animal's session compiled data
else 
    load([dir filesep data(m).expID{end} filesep 'compiled analyzed data_' hva '.mat']); % load last animal's session compiled data 
end 
% for exp = 1: size(analyzed_comp,1) % find each experiment's number of contrasts
%     if gratings 
%         cont(exp) = length(analyzed_comp(exp,1).contrasts);
%     else 
%         cont(exp) = length(analyzed_comp(exp,1).cont1);
%     end
% end 
% nrContrasts = max(cont);

%%
nrContrasts = length(data(1).idxb);
temp = cell(1,size(data(1).avg_exp,2));
for m =1: length(data) 
    for i = 1: size(data(m).avg_exp,2) % locations 
        temp{1,i} = cat(3, temp{1,i}, data(m).avg_exp{1,i}); % concatenate over animals 
 
    end
end 

avg_animal = cellfun(@(x) nanmean(x, 3), temp, 'UniformOutput', false);
std_animal =  cellfun(@(x) nanstd(x, 0, 3), temp, 'UniformOutput', false);

for j = 1: size(avg_animal, 2) % locations
    for i = 1:nrContrasts
        [peak_value(i,j), peak_index] = getPeakPostStim(avg_animal{j}(i,:), 15, 1, 2);
        std_peak(i,j) = std_animal{j}(i, peak_index);
    end
end 
% end 

f = figure;
subplot(1,2,1)
errorbar(data(1).idxb*100, peak_value(:,1)*100, std_peak(:,1)*100, 'ko-', 'LineWidth', 1.5);
hold on; 

xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('Binocular'); axis square; axis padded; box off

subplot(1,2,2)
errorbar(data(1).idxm*100, peak_value(:,2)*100, std_peak(:,2)*100, 'ko-', 'LineWidth',1.5);

xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('Monocular'); axis square; axis padded; box off
sgtitle(['average response activity across animals (N = ' num2str(length(Animal)) ')']);

% yLim = [-0.01 0.035];
% yLim = [min(min(peak_value))-max(max(std_peak)) max(max(peak_value))+max(max(std_peak))]; % Define the y-axis limit
% ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
% set(ax,'ylim',yLim);

%% data is saved nrContrasts x frames x sessions across all experiments 
dir = [fPath 'stats'];
if ~isfolder(dir)
    mkdir (dir)
end
dir = [dir filesep 'behavior'];
if ~isfolder(dir)
    mkdir (dir)
end

cd(dir)
if isempty(hva)
    savefig(f, ['meandff_over_contrasts_animals_' type '.fig'])
    saveas(f, ['meandff_over_contrasts_animals_' type '.png'])
    save(['compiled ROI traces_animals_' type '.mat'], 'temp', 'Animal', 'data')
else
    savefig(f, ['meandff_over_contrasts_animals_' type '_' hva '.fig'])
    saveas(f, ['meandff_over_contrasts_animals_' type '_' hva '.png'])
    save(['compiled ROI traces_animals_' type '_' hva '.mat'], 'temp', 'Animal', 'data')
end
%% 
figure; 
subplot 121
[traceb,stdb] = getAvgTraceAcrossStimLocs(temp(1), data(1).idxb, 30, 0);
title('Binocular');
subplot 122
[tracem,stdm] = getAvgTraceAcrossStimLocs(temp(2), data(1).idxm, 30, 0);
title('Monocular');
sgtitle ('Average Response Across Animals')

yLim = [-1.5 4];
ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
set(ax,'ylim',yLim);
if isempty(hva)
    savefig(gcf, ['meandffTrace_animals_' type '.fig'])
    saveas(gcf, ['meandffTrace_animals_' type '.png'])
else
    savefig(gcf, ['meandffTrace_animals_' type '_' hva '.fig'])
    saveas(gcf, ['meandffTrace_animals_' type '_' hva '.png'])
end

%%
% figure;
% for i = 1:size(trace,1)
%     [peak_cont(i), peak_index] = getPeakPostStim(trace(i,:), 15, 1, 2);
%     std_peak_cont(i) = std(i, peak_index);
% end
% 
% errorbar(cont, peak_cont*100, std_peak_cont*100, 'o');
% xlabel('Contrasts (%)')
% ylabel('\Deltaf/f (%)');
% title('avg response across stimulus locations'); axis square; axis padded; box off
% ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
% % set(ax,'ylim',[-0.2 2]);
% 
% savefig(gcf, 'meandff_over_contrasts_animals_allLocations.fig')
% saveas(gcf, 'meandff_over_contrasts_animals_allLocations.png')