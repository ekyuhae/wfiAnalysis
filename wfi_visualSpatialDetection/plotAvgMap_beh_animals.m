%%plotAvgMap_beh_days across animals and do hit-miss
%% load all hit data
[binoc_array_1, monoc_array_1] = plotAvgMap_beh_days('M230220_1','hit','RL_AM',0);
[binoc_array_2, monoc_array_2] = plotAvgMap_beh_days('M230209_1','hit','RL_AM',0);


addpath('Y:\\haider\Data\analyzedData\EKK\WFI\M230220_1\retinotopy\04-Apr-2023_1');
load('Y:\\haider\Data\analyzedData\EKK\WFI\M230220_1\retinotopy\04-Apr-2023_1\plotPhaseMap.mat')
plotPhaseMap1 = plotPhaseMap;

addpath('Y:\\haider\Data\analyzedData\EKK\WFI\M230209_1\retinotopy\15-Mar-2023')
load('Y:\\haider\Data\analyzedData\EKK\WFI\M230209_1\retinotopy\15-Mar-2023\plotPhaseMap.mat')
plotPhaseMap2 = plotPhaseMap;

load('Y:\\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK\allenDorsalMapSM.mat')
%%

%DO FOR BOTH ANIMALS - image align to allen and iterate through contrasts
moved1 = cell(1,6);
[combined, moved, fixed,cord] = imalign_AR(plotPhaseMap1, dorsalMaps.areaMap); %for selection of contour onto allen
movedA1 = moved;
for c = 1:6 %iterate through contrasts
    [moved,combined] = imalign_s(binoc_array_1{c}, fixed,cord.moved,cord.fixed); %for avgImage onto already selected align points
    moved1{c} = moved;
end

moved2 = cell(1,6);
[combined, moved, fixed,cord] = imalign_AR(plotPhaseMap2, dorsalMaps.areaMap); %for selection of contour onto allen
movedA2 = moved;
for c = 1:6 % iterate through contrasts
    [moved,combined] = imalign_s(binoc_array_2{c}, fixed,cord.moved,cord.fixed); %for avgImage onto already selected align points
    moved2{c} = moved;
end

%%
close all
clearvars -except dorsalMaps binoc_array_1 monoc_array_1 binoc_array_2 monoc_array_2 plotPhaseMap1 plotPhaseMap2 moved1 moved2
%%
f1 = figure;
for c = 1:6
    movedT = (moved1{c} + moved2{c})./2;
    %f = figure; imagesc(movedT);
    fans = movedT(250:550,1:350);
    %f2 = figure; imagesc(fans)
    subplot(3,2,c, 'Parent', f1); 
    imagesc(fans);
    hold on
end

%use imalign to align contours
%% average hit trials across animals binoc
for i = 1:length(binoc_array_1)
    binoc_array_1{i} = cat(3,binoc_array_1{i}, binoc_array_2{i});
end
binoc_array = binoc_array_1;

cmin = []; cmax = [];
for i = 1:length(binoc_array)
    binoc_array{i} = mean(binoc_array{i},3);
    cmin = [cmin min(min(binoc_array{i}))];
    cmax = [cmax max(max(binoc_array{i}))];
end
cmin = min(cmin); cmax = max(cmax);
clim([cmin cmax]); hold on
for i = 1:length(bcont) %plot
    binoc_array{i} = mean(binoc_array{i},3);
    nexttile
    imagesc(binoc_array{i}); colormap viridis; freezeColors; colorbar; clim([cmin cmax])
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    if i > 1
        title(sprintf('%d%%', bcont(i)*100));
    else
        title(typeo)
    end
end
sgtitle(['avg activity map across animals binoc- RL AM hit v FA trials'])
set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(m, ['Y:\haider\Data\analyzedData\EKK\WFI\avgActivityMap_binoc_hitTrials.fig'])
binoc_array = binoc_array_hit;

%% average hit trials across animals monoc
for i = 1:length(monoc_array_1)
    monoc_array_1{i} = cat(3,monoc_array_1{i}, monoc_array_2{i});
end
monoc_array = monoc_array_1;

cmin = []; cmax = [];
for i = 1:length(monoc_array)
    monoc_array{i} = mean(monoc_array{i},3);
    cmin = [cmin min(min(monoc_array{i}))];
    cmax = [cmax max(max(monoc_array{i}))];
end
cmin = min(cmin); cmax = max(cmax);
clim([cmin cmax]); hold on
for i = 1: length(mcont)
    nexttile
    imagesc(monoc_array{i}); colormap viridis; freezeColors; colorbar; clim([cmin cmax]);
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    if i > 1
        title(sprintf('%d%%', mcont(i)*100));
    else
        title(typeo)
    end
end
sgtitle(['avg activity map across days monoc- RL AM hit v FA trials'])
set(f, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(f, ['Y:\haider\Data\analyzedData\EKK\WFI\avgActivityMap_monoc_hitTrials.fig'])
monoc_array = monoc_array_hit;

%% load all miss data
[binoc_array_1, monoc_array_1] = plotAvgMap_beh_days('M230220_1','miss','RL_AM',0);
[binoc_array_2, monoc_array_2] = plotAvgMap_beh_days('M230209_1','miss','RL_AM',0);

%% average miss trials across animals binoc
for i = 1:length(binoc_array_1)
    binoc_array_1{i} = cat(3,binoc_array_1{i}, binoc_array_2{i});
end
binoc_array = binoc_array_1;

cmin = []; cmax = [];
for i = 1:length(binoc_array)
    binoc_array{i} = mean(binoc_array{i},3);
    cmin = [cmin min(min(binoc_array{i}))];
    cmax = [cmax max(max(binoc_array{i}))];
end
cmin = min(cmin); cmax = max(cmax);
clim([cmin cmax]); hold on
for i = 1:length(bcont) %plot
    binoc_array{i} = mean(binoc_array{i},3);
    nexttile
    imagesc(binoc_array{i}); colormap viridis; freezeColors; colorbar; clim([cmin cmax])
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    if i > 1
        title(sprintf('%d%%', bcont(i)*100));
    else
        title(typeo)
    end
end
sgtitle(['avg activity map across animals binoc- RL AM miss v CR trials'])
set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(m, ['Y:\haider\Data\analyzedData\EKK\WFI\avgActivityMap_binoc_missTrials.fig'])
binoc_array = binoc_array_miss;

%% average miss trials across animals monoc
for i = 1:length(monoc_array_1)
    monoc_array_1{i} = cat(3,monoc_array_1{i}, monoc_array_2{i});
end
monoc_array = monoc_array_1;

cmin = []; cmax = [];
for i = 1:length(monoc_array)
    monoc_array{i} = mean(monoc_array{i},3);
    cmin = [cmin min(min(monoc_array{i}))];
    cmax = [cmax max(max(monoc_array{i}))];
end
cmin = min(cmin); cmax = max(cmax);
clim([cmin cmax]); hold on
for i = 1: length(mcont)
    nexttile
    imagesc(monoc_array{i}); colormap viridis; freezeColors; colorbar; clim([cmin cmax]);
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    if i > 1
        title(sprintf('%d%%', mcont(i)*100));
    else
        title(typeo)
    end
end
sgtitle(['avg activity map across days monoc- RL AM miss v CR trials'])
set(f, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(f, ['Y:\haider\Data\analyzedData\EKK\WFI\avgActivityMap_monoc_missTrials.fig'])
monoc_array = monoc_array_miss;

%% Hits - Miss binoc

for i = 1:length(binoc_array_hit)
    binoc_array_hit{i} = binoc_array_hit{i} - binoc_array_miss{i};
end
binoc_array = binoc_array_hit;

cmin = []; cmax = [];
for i = 1:length(binoc_array)
    cmin = [cmin min(min(binoc_array{i}))];
    cmax = [cmax max(max(binoc_array{i}))];
end
cmin = min(cmin); cmax = max(cmax);
clim([cmin cmax]); hold on
for i = 1:length(bcont) %plot
    nexttile
    imagesc(binoc_array{i}); colormap viridis; freezeColors; colorbar; clim([cmin cmax])
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    if i > 1
        title(sprintf('%d%%', bcont(i)*100));
    else
        title(typeo)
    end
end
sgtitle(['avg activity map across animals binoc hits - miss'])
set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(m, ['Y:\haider\Data\analyzedData\EKK\WFI\avgActivityMap_binoc_hitsminusmiss.fig'])

%% Hits - miss monoc
for i = 1:length(monoc_array_hit)
    monoc_array_hit{i} = monoc_array_hit{i}- monoc_array_miss{i};
end
monoc_array = monoc_array_hit;

cmin = []; cmax = [];
for i = 1:length(monoc_array)
    cmin = [cmin min(min(monoc_array{i}))];
    cmax = [cmax max(max(monoc_array{i}))];
end
cmin = min(cmin); cmax = max(cmax);
clim([cmin cmax]); hold on
for i = 1: length(mcont)
    nexttile
    imagesc(monoc_array{i}); colormap viridis; freezeColors; colorbar; clim([cmin cmax]);
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    if i > 1
        title(sprintf('%d%%', mcont(i)*100));
    else
        title(typeo)
    end
end
sgtitle(['avg activity map across days monoc hits - miss'])
set(f, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
savefig(f, ['Y:\haider\Data\analyzedData\EKK\WFI\avgActivityMap_monoc_hitsminusmiss.fig'])
monoc_array = monoc_array_miss;


