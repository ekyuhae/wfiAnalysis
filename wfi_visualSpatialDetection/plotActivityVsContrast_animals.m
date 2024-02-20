%% load individual session datasets that were analyzed from lv2 and compile them 
% modified version of stat_new by AR
% EK 23 

close all; clear all; clc;
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
% addpath(genpath('\\ad.gatech.edu\bme\labs\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\';
Animal = {'M221118_1', 'M230130_2', 'M230131_1'};
% expID = {'03-Feb-2023', '08-Feb-2023_2','18-Feb-2023', '09-Mar-2023'};
expID = {'09-Feb-2023', '16-Feb-2023', '18-Feb-2023', '26-Feb-2023', '27-Feb-2023', '02-Mar-2023','02-Mar-2023_1'};
gratings = true; % false if barmapping
% cont = [100 75 50 25 5 0]; % barmapping contrast (shortcut)
cont = [0,0.02,0.05,0.1,0.33,0.5,0.65,0.8,1]*100; % grating cont
%%
for m = 1: length(Animal)
    if gratings
        dir = [fPath Animal{m} filesep 'gratings'];
        cd(dir)
        data(m) = load('compiled ROI traces_days.mat');

    else 
          dir = [fPath Animal{m} filesep 'barmapping'];
        cd(dir)
        data(m) = load('compiled ROI traces_days.mat');
    end 
     
end 
load([dir filesep data(m).expID{end} filesep 'compiled analyzed data.mat']); % load last animal's session compiled data 

% for exp = 1: size(analyzed_comp,1) % find each experiment's number of contrasts
%     if gratings 
%         cont(exp) = length(analyzed_comp(exp,1).contrasts);
%     else 
%         cont(exp) = length(analyzed_comp(exp,1).cont1);
%     end
% end 
% nrContrasts = max(cont);

%%
nrContrasts = size(data(1).avg_exp{1},1);
temp = cell(1,size(data(1).avg_exp,2));
for m =1: length(data) 
    for i = 1: size(data(m).avg_exp,2) % locations 
        temp{1,i} = cat(3, temp{1,i}, data(m).avg_exp{1,i}); % concatenate over sessions
 
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
if gratings 
errorbar(analyzed_comp(end,1).contrasts*100, peak_value(:,1)*100, std_peak(:,1)*100, 'ko' ,'LineWidth', 1.5);
else 
    errorbar(cont, peak_value(:,5), std_peak(:,5), 'ko');
end 
xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('Binocular V1'); axis square; axis padded; box off

subplot(1,2,2)
if gratings 
errorbar(analyzed_comp(end,1).contrasts*100, peak_value(:,2)*100, std_peak(:,2)*100, 'ko', 'LineWidth', 1.5);
else 
    errorbar(cont, peak_value(:,12), std_peak(:,12), 'o');
end
xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('Monocular V1'); axis square; axis padded; box off
sgtitle(['average response activity across animals (N = ' num2str(length(Animal)) ')']);

% yLim = [-0.01 0.035];
% % yLim = [min(min(peak_value))-max(max(std_peak)) max(max(peak_value))+max(max(std_peak))]; % Define the y-axis limit
% ax = findobj(f,'type','axes'); % Find all axes handles in the figure
% set(ax,'ylim',yLim);

%% data is saved nrContrasts x frames x sessions across all experiments 
dir = [fPath 'stats'];
if ~isfolder(dir)
    mkdir (dir)
end
if gratings
    dir = [dir filesep 'gratings'];
    if ~isfolder(dir)
        mkdir (dir)
    end
else
    dir = [dir filesep 'barmapping'];
    if ~isfolder(dir)
        mkdir (dir)
    end
end
cd(dir)
savefig(f, 'meandff_over_contrasts_animals.fig')
saveas(f, 'meandff_over_contrasts_animals.png')
save('compiled ROI traces_animals.mat', 'temp', 'Animal')

%% 

[trace,std] = getAvgTraceAcrossSess(temp, cont/100, 30,1);

figure;
for i = 1:size(trace,1)
    [peak_cont(i), peak_index] = getPeakPostStim(trace(i,:), 15, 1, 2);
    std_peak_cont(i) = std(i, peak_index);
end

errorbar(cont, peak_cont*100, std_peak_cont*100, 'ko', 'LineWidth', 1.5);
xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('avg response across stimulus locations'); axis square; axis padded; box off
ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
% set(ax,'ylim',[-0.2 2]);

savefig(gcf, 'meandff_over_contrasts_animals_allLocations.fig')
saveas(gcf, 'meandff_over_contrasts_animals_allLocations.png')