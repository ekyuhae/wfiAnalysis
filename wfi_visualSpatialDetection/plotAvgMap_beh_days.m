function [binoc_array, monoc_array] = plotAvgMap_beh_days(Animal, type, hva, reactionTriggered)

% average activity maps per contrast across all days in one animal
% INPUTS: [1] Animal = 'M230220_1'; [2] type = 'hit'; [3] hva = 'RL_AM'; [4] reactionTriggered = 0;
% OUTPUTS: 2 figures (one for binoc region ROI and one for monoc region ROI), w spatiotemporal response for 6 contrasts
% similar to: plotAvgMap_beh.m and plotAvgMap_beh_animals.m

% Used in: plotAvgMap_beh_animals.m

% NA 6/22/2023, 6/28

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
aPath = ['Y:\haider\Data\analyzedData\EKK\WFI' filesep Animal filesep 'behavior'];

if strcmp(type, 'hit')
    % for hit trials, false alarms are used for the 0% (licking)
    typeo = 'FA';
elseif strcmp(type, 'miss')
    % for miss trials, correct rejects are used for the 0% (no licking)
    typeo = 'CR';
end

% generate list of all days for the animal
list = dir(aPath); explist = []; 
for i = 1:height(list)
    if double(list(i,1).name(1)) >= 48 && double(list(i,1).name(1) <=57) %only the dates
        explist = [explist; {list(i,1).name}];
    end
end

explist ={'14-Aug-2023' ;'15-Aug-2023';'17-Aug-2023'};

avgStimRA.binoc = cell(length(explist),6); avgStimRA.monoc = cell(length(explist),6); %initialize data structs
for d = 1:length(explist) % iterate through days
    dayDir = dir([aPath filesep explist{d,:}]);
    % generate list of all sessions for that day
    sessList = [];
    for ii = 1:height(dayDir)
        if double(dayDir(ii).name(1) >= 48 && double(dayDir(ii).name(1)<=57)) && length(dayDir(ii).name) <=2
            sessList = [sessList; dayDir(ii).name];
        end
    end

    for v = 1:2 % iterate between visual areas (binoc, monoc)
        if v == 1
            mb = 'binoc';
        elseif v == 2
            mb = 'monoc';
        end

        for s = 1:length(sessList) % iterate through sessions
            try % skip things that don't have data
                sPath = [aPath filesep explist{d,:} filesep sessList(s,:)];
                load([sPath filesep 'RL_AM' filesep 'analyzed_active_' num2str(s) '_' hva '_' type '.mat'])
                for ss = 2:length(avgStimResponse.(mb)) % iterate through contrasts but not 0%
                    if isempty(avgStimResponse.(mb){1,ss})
                        avgStimRA.(mb){d,ss} = [];
                    else
                        avgStimRA.(mb){d,ss} = cat(3,avgStimRA.(mb){d,ss}, avgStimResponse.(mb){1,ss});
                    end
                end
                % add in 0% contrasts data (FA if hits, CR if miss)
                load([sPath filesep 'RL_AM' filesep 'analyzed_active_' num2str(s) '_' hva '_' typeo '.mat'])
                if isempty(avgStimResponse.(mb){1,1})
                    avgStimRA.(mb){d,1} = [];
                else
                    avgStimRA.(mb){d,1} = cat(3, avgStimRA.(mb){d,1}, avgStimResponse.(mb){1,1});
                end
            end
        end

        for ss = 1:width(avgStimRA.(mb))  % average across sessions
            avgStimRA.(mb){d,ss} = mean(avgStimRA.(mb){d,ss},3);
        end
    end
end

clearvars -except avgStimRA contrasts explist Animal type typeo

% then iterate the resize and align across days - COPIED from plotAvgMap_beh.m
load(['Y:\haider\Data\analyzedData\EKK\WFI\' Animal filesep 'behavior' filesep explist{1,:} '\1\contours.mat'])
contours = imread(['Y:\haider\Data\analyzedData\EKK\WFI\' Animal filesep 'behavior' filesep explist{1,:}  '\1\overlaid_contourMap.png']);
contours = rgb2ind(contours, 256);
idx = find(contrasts == 0, 1, 'last');

% plotting
for v = 1:2 % iterate between binoc and monoc
    if v == 1
        mb = 'binoc';
        cont = contrasts(1:idx-1);
    elseif v == 2
        mb = 'monoc';
        cont = contrasts(idx:end);
    end

    for d = 1: height(avgStimRA.(mb)) % iterate through days
        for c = 1: width(avgStimRA.(mb))% iterate through contrasts
            if ~isempty(avgStimRA.(mb){d,c})
                data{d,c} = avgStimRA.(mb){d,c};
            else
                data{d,c}= NaN;
            end
        end
    end
    m = figure;
    tiledlayout(2, (1/2)*length(cont),'TileSpacing', 'compact', 'Padding', 'compact'); % show avg stim response activity for each white bar location
    data_array = cell(1,6);
    for d = 1:height(avgStimRA.(mb)) % iterate through days
        for i = 1:length(cont) % iterate through contrasts
            data_resized{d,i} = imresize(data{d,i}, [size(contours,1), size(contours,2)]);
            load(['Y:\haider\Data\analyzedData\EKK\WFI\' Animal filesep 'behavior' filesep explist{d,:} '\1\coordinates.mat'])
            [data_resized{d,i},~] = imalign_s(data_resized{d,i}, contours, cord.moved,cord.fixed);
            if ~isnan(sum(sum(data_resized{d,i}))) % avoid sessions w/out data
                data_array{i} = cat(3,data_array{i}, data_resized{d,i});
            end
        end
    end
    cmin = []; cmax = [];
    for i = 1:length(data_array)
        data_array{i} = mean(data_array{i},3);
        cmin = [cmin min(min(data_array{i}))];
        cmax = [cmax max(max(data_array{i}))];
    end
    cmin = min(cmin); cmax = max(cmax);
    clim([cmin cmax]); hold on
    for i = 1:length(cont) % plot
        data_array{i} = mean(data_array{i},3);
        nexttile
        imagesc(data_array{i}); colormap viridis; freezeColors; colorbar; clim([cmin cmax])
        hold on
        contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
        contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
        axis image
        if i > 1
            title(sprintf('%d%%', cont(i)*100));
        else
            title(typeo)
        end
    end
    sgtitle(['avg activity map across days ' mb ' - RL AM ' type ' v ' typeo ' trials'])
    set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
    savefig(m, ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior\avgActivityMap_' mb '_' type 'Trials.fig'])
    if v == 1
        binoc_array = data_array;
    elseif v ==2
        monoc_array = data_array;
    end
end
end