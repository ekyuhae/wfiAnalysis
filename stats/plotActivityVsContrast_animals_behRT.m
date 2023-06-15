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

type = 'hit';
hva = 'RL_AM';
RT = '_RT';
%% hits
for m = 1: length(Animal)
    dir = [fPath Animal{m} filesep 'behavior'];
    if isempty(hva)
        stimTrig(m) = load([dir filesep 'V1_compiled_days' filesep 'compiled ROI traces_days_' type '.mat']);
        lickTrig(m) = load([dir filesep 'V1_RT_compiled_days' filesep 'compiled ROI traces_days_' type '.mat']);
    else
        stimTrig(m) = load([dir filesep hva '_compiled_days' filesep 'compiled ROI traces_days_' type '_' hva '.mat']);
        lickTrig(m) = load([dir filesep hva RT '_compiled_days' filesep 'compiled ROI traces_days_' type '_' hva '.mat']);
    end
end

nrContrasts = length(stimTrig(1).idxb);
for m = 1:length(Animal)
    mags(m).binoc = getPeakEvokedAmp(lickTrig(m).avg_exp{1}, stimTrig(m).avg_exp{1}, 15, 2);
    mags(m).monoc = getPeakEvokedAmp(lickTrig(m).avg_exp{2}, stimTrig(m).avg_exp{2}, 15, 2);
end

mags_all_b = [mags(1).binoc mags(2).binoc];
mags_all_m = [mags(1).monoc mags(2).monoc];

peak_value(:,1) = nanmean(mags_all_b,2); 
peak_value(:,2) = nanmean(mags_all_m,2); 
std_peak(:,1) = nanstd(mags_all_b,0,2);
std_peak(:,2) = nanstd(mags_all_m,0,2);
%% misses

for m = 1: length(Animal)
    dir = [fPath Animal{m} filesep 'behavior'];
    if isempty(hva)
        stimTrig_miss(m) = load([dir filesep 'V1_compiled_days' filesep 'compiled ROI traces_days_miss.mat']);
    else
        stimTrig_miss(m) = load([dir filesep hva '_compiled_days' filesep 'compiled ROI traces_days_miss_' hva '.mat']);
    end
end

for m = 1:length(Animal)
    ST = [];
    ST = stimTrig_miss(m).avg_exp{1};
    for iSess = 1: size(stimTrig_miss(m).avg_exp{1},3) %loop through sessions
        for iCont = 1: size (stimTrig_miss(m).avg_exp{1},1) %loop through contrasts
            [peak_val, ~] = getPeakPostStim(ST(iCont, :, iSess), 15, 1, 2, []);
            magnitude(iCont, iSess) = peak_val;
        end
    end
    mags_miss(m).binoc = magnitude;

    ST = stimTrig_miss(m).avg_exp{2};
    for iSess = 1: size(stimTrig_miss(m).avg_exp{2},3) %loop through sessions
        for iCont = 1: size (stimTrig_miss(m).avg_exp{2},1) %loop through contrasts
            [peak_val, ~] = getPeakPostStim(ST(iCont, :, iSess), 15, 1, 2, []);
            magnitude(iCont, iSess) = peak_val;
        end
    end
    mags_miss(m).monoc = magnitude;
    clearvars magnitude
end

mags_all_miss_b = [mags_miss(1).binoc mags_miss(2).binoc];
mags_all_miss_m = [mags_miss(1).monoc mags_miss(2).monoc];

peak_value_miss(:,1) = nanmean(mags_all_miss_b,2); 
peak_value_miss(:,2) = nanmean(mags_all_miss_m,2); 
std_peak_miss(:,1) = nanstd(mags_all_miss_b,0,2);
std_peak_miss(:,2) = nanstd(mags_all_miss_m,0,2);

%% figure with miss

color = [0.77,0.44,0.44];
colorM = [128 128 128]./255;

f = figure;
subplot(1,2,1)
errorbar(stimTrig_miss(1).idxb*100, peak_value_miss(:,1)*100, std_peak_miss(:,1)*100, 'o-', 'color', colorM, 'LineWidth', 1.5, 'CapSize',0,'MarkerFaceColor',colorM);
hold on; 
errorbar(stimTrig(1).idxb*100, peak_value(:,1)*100, std_peak(:,1)*100, 'o-', 'color', color, 'LineWidth', 1.5, 'CapSize',0,'MarkerFaceColor',color); 

xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
if isempty(hva)
    title('Binocular V1');
else
    title (hva(1:2))
end
legend('Miss', 'Hit', 'Location', 'northwest'); legend boxoff
axis square; axis padded; box off

subplot(1,2,2)
errorbar(stimTrig_miss(1).idxm*100, peak_value_miss(:,2)*100, std_peak_miss(:,2)*100, 'o-', 'color', colorM, 'LineWidth', 1.5, 'CapSize',0,'MarkerFaceColor',colorM);
hold on;
errorbar(stimTrig(1).idxm*100, peak_value(:,2)*100, std_peak(:,2)*100, 'o-', 'color', color, 'LineWidth', 1.5, 'CapSize',0,'MarkerFaceColor',color); 

xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
if isempty(hva)
    title('Monocular V1');
else
    title (hva(4:5))
end
axis square; axis padded; box off
sgtitle(['Average response activity across animals (N = ' num2str(length(Animal)) ')']);
legend('Miss', 'Hit', 'Location', 'northwest');  legend boxoff

yLim = [-1 3.5];
% yLim = [min(min(peak_value))-max(max(std_peak)) max(max(peak_value))+max(max(std_peak))]; % Define the y-axis limit
ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
set(ax,'ylim',yLim);
%%
dir = [fPath 'stats' filesep 'behavior'];
cd(dir)
mags_all.hits.binoc = mags_all_b; 
mags_all.hits.monoc = mags_all_m;
mags_all.miss.binoc = mags_all_miss_b;
mags_all.miss.monoc = mags_all_miss_m;

if isempty(hva)
    savefig(f, ['dFF_hit_miss.fig'])
    saveas(f, ['dFF_hit_miss.png'])
   
else
    savefig(f, ['dFF_hit_miss_' hva '.fig'])
    saveas(f, ['dFF_hit_miss_' hva '.png'])
   
end

%% comparison of peak sensory responses between detection hit and miss trials 

x_labels = {'Miss', 'Hit'};
x = 1:length(x_labels);

figure; 
subplot 121
k = linspace(0.9,0.2, nrContrasts);
for i = 1:nrContrasts
    c = [k(i), k(i), k(i)];
    errorbar(x, [peak_value_miss(i,1)*100 peak_value(i,1)*100], [std_peak_miss(i,1)*100 std_peak(i,1)*100], 'o-', 'CapSize',0,'color',c,'MarkerFaceColor',c); 
    hold on
end
set(gca, 'XTick', x, 'XTickLabel', x_labels);
axis padded; axis square; box off
set(findobj(gcf, 'type', 'axes'), 'xlim', [0 3], 'ylim', [-1 3])
legend('0%' , '1%', '2%', '5%', '10%', '18%', 'Location','northwest'); legend boxoff
ylabel('Mean Peak Value (%)');
if isempty(hva)
    title('Binocular V1');
else
    title (hva(1:2))
end

for i = 1:nrContrasts
    [~, Pvals(i,1)] = ttest(mags_all_b(i,:)', mags_all_miss_b(i,1:size(mags_all_b,2))');
end

% Add p-values as text annotations
p_value_threshold = 0.05;
for i = 1:nrContrasts
    if ~isnan(Pvals(i,1)) && Pvals(i,1) < p_value_threshold
        text(x(2) + 0.1, peak_value(i,1)*100, "*", 'FontSize', 12); hold on;
    end
end

subplot 122
for i = 1:nrContrasts
    c = [k(i), k(i), k(i)];
    errorbar(x, [peak_value_miss(i,2)*100 peak_value(i,2)*100], [std_peak_miss(i,2)*100 std_peak(i,2)*100], 'o-', 'CapSize',0,'color',c,'MarkerFaceColor',c); 
    hold on
end
set(gca, 'XTick', x, 'XTickLabel', x_labels);
axis padded; axis square; box off
set(findobj(gcf, 'type', 'axes'), 'xlim', [0 3],  'ylim', [-1 3])
legend('0%' , '2%', '5%', '10%', '18%', '33%','Location','northwest'); legend boxoff
ylabel('Mean Peak Value (%)');
if isempty(hva)
    title('Monocular V1');
else
    title (hva(4:5))
end

for i = 1:nrContrasts
    [~, Pvals(i,2)] = ttest(mags_all_m(i,:)', mags_all_miss_m(i,1:size(mags_all_m,2))');
end

p_value_threshold = 0.05;
for i = 1:nrContrasts
    if ~isnan(Pvals(i,2)) && Pvals(i,2) < p_value_threshold
        text(x(2) + 0.1, peak_value(i,2)*100, "*", 'FontSize', 12); hold on;
    end
end
%% Perform paired t-test between the groups for each contrast level

if isempty(hva)
    savefig(gcf, 'CompareHitVsMiss.fig');
    saveas(gcf, 'CompareHitVsMiss.png');
    save('respAmpComparison.mat', 'mags_all', 'Pvals');
else
    savefig(gcf, ['CompareHitVsMiss_' hva '.fig']);
    saveas(gcf, ['CompareHitVsMiss_' hva '.png']);
    save(['respAmpComparison' hva '.mat'], 'mags_all', 'Pvals');
end