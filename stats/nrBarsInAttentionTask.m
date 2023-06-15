close all; clear all; clc;
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

fPath = 'Y:\haider\Data\Behavior\expInfo';
Animal = {'M201210_1', 'M220502_1' 'M220513_1' };
expID = {'2021-03-10', '2022-06-30' '2022-06-23' };

%%
for i = 1:length(Animal)
    cPath = [fPath filesep Animal{i} filesep expID{i}];
       
         % get a list of all subfolders
         sessions = dir(cPath);
         sessions = sessions([sessions(:).isdir]); % keep only directories
         sessions = sessions(3:end); % remove '.' and '..'
         
         % count the number of subfolders whose names contain only numbers
%          count = 0;
%          for i = 1:length(sessions)
%              if regexp(sessions(i).name, '^\d+$')
%                  count = count + 1;
%              end
%          end
         nrSessions(i) = length(sessions);
         
         % compile analyzed mat. file from each session (column) and compile
         % them with different experiments
         for sessNr = 1:nrSessions(i)
             p = [cPath filesep num2str(sessNr)];
             cd(p)
             cFile =[expID{i} '_' num2str(sessNr) '_' Animal{i} '_Block.mat'];
             if exist(cFile)
                blocks(i, sessNr)= load([expID{i} '_' num2str(sessNr) '_' Animal{i} '_Block.mat']);
             end 
         end 
end 

%%
for i = 1:size(blocks,1)
    
    
    for s = 1: nrSessions(i)
        
        if ~isempty(blocks(i,s).block)
            flag = strcmp(blocks(i,s).block.expType, 'BarMapping');
            if ~flag               
                for k = 1: length(blocks(i,s).block.trial)
                    b(k) = blocks(i,s).block.trial(k).condition.targetAzimuth <= 45;
                    m(k) = blocks(i,s).block.trial(k).condition.targetAzimuth >= 45; 
                end  
                locs{i}(s,b==1) = 1; % idx for binoc block 
                locs{i}(s,m==1) = 2; 
                b = []; m = [];
            else 
                locs{i}(s,:) = NaN;
            end 
        else 
            locs{i}(s,:) = NaN;
        end
    end
    
end

%% 
for i = 1:size(blocks,1)
     for s = 1: nrSessions(i)
        temp= [];
        monoc = find(locs{i}(s,:) == 2);
        if ~isempty(monoc) && strcmp(blocks(i,s).block.endStatus, 'completed')
            range = [blocks(i,s).block.trial.trialStartedTime; blocks(i,s).block.trial.quiescentEpochTime];
            for m = monoc
                idx = find(blocks(i,s).block.barChangeTime >= range(1,m) & blocks(i,s).block.barChangeTime <= range(2,m));
                count = length(idx);
                temp = [temp count];
            end
            nrBar_monoc{i, s} = temp;
        end
    end
end

for i = 1:size(blocks,1)
    for s = 1: nrSessions(i)
        temp= [];
        binoc = find(locs{i}(s,:) == 1);
        if ~isempty(binoc) && strcmp(blocks(i,s).block.endStatus, 'completed')
            range = [blocks(i,s).block.trial.trialStartedTime; blocks(i,s).block.trial.quiescentEpochTime];
            for b = binoc
                idx = find(blocks(i,s).block.barChangeTime >= range(1,b) & blocks(i,s).block.barChangeTime <= range(2,b));
                count = length(idx);
                temp = [temp count];
            end
            nrBar_binoc{i, s} = temp;
        end
    end
end

%% 
% 22 blocks for binoc/monoc grating, 
%15 trials for binoc, 30 trials for monoc 

nrBar_binoc = reshape(nrBar_binoc', 1,[]);
nrBarAll_binoc = [];
for i =1:length(nrBar_binoc)
nrBarAll_binoc = [nrBarAll_binoc nrBar_binoc{i}];
end 

nrBar_monoc = reshape(nrBar_monoc', 1,[]);
nrBarAll_monoc = [];
for i =1:length(nrBar_monoc)
nrBarAll_monoc = [nrBarAll_monoc nrBar_monoc{i}];
end 

avgNr_b= nanmean(nrBarAll_binoc);
avgNr_m = nanmean(nrBarAll_monoc);
std_b = nanstd(nrBarAll_binoc);
std_m = nanstd(nrBarAll_monoc);
max_b = max(nrBarAll_binoc); max_m = max(nrBarAll_monoc);
min_b = min(nrBarAll_binoc); min_b = min(nrBarAll_monoc);

figure; errorbar(0, avgNr_b, std_b, 'ok'); hold on; errorbar(70, avgNr_m, std_m, 'ok')

% boxplot({nrBarAll_monoc, nrBarAll_binoc}); hold on
b = cell(length(nrBarAll_binoc),1); m = cell(length(nrBarAll_monoc),1);
b(:) = {'binoc'}; m(:)={'monoc'};
group = [b;m];
figure; 
boxx = boxchart(categorical(group), [nrBarAll_binoc'; nrBarAll_monoc']); 
hold on;
plot([1,2], [mean(nrBarAll_binoc) mean(nrBarAll_monoc)], 'o', 'Color', [0.5 0.5 1]);
boxx.JitterOutliers = 'on';
boxx.MarkerStyle = '.';
ylim([0 300]); 
box off; axis square
% xlabel ('degree'); 
ylabel('Nr of Bars within each block')
