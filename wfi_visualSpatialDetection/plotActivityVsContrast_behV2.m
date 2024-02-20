close all; clear; clc


addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK\project'))
addpath(genpath("Y:\haider\Code\behaviorAnalysis"));

Animals = {'M230717_1'}; %cell array
fPath = 'Y:\haider\Data\analyzedData\EKK\WFI';
% expID = {'14-Sep-2023',  '15-Sep-2023' '16-Sep-2023', '17-Sep-2023', '18-Sep-2023_1', '19-Sep-2023', '20-Sep-2023', '21-Sep-2023','22-Sep-2023', '26-Sep-2023', '27-Sep-2023', '29-Sep-2023', '01-Oct-2023', '02-Oct-2023', '04-Oct-2023','04-Oct-2023', '06-Oct-2023', '07-Oct-2023_1', '09-Oct-2023', '12-Oct-2023','13-Oct-2023', '15-Oct-2023', '15-Oct-2023_1','17-Oct-2023','19-Oct-2023', '20-Oct-2023'};
%expID = {'14-Sep-2023',  '15-Sep-2023' '16-Sep-2023', '17-Sep-2023', '18-Sep-2023_1', '19-Sep-2023', '20-Sep-2023', '21-Sep-2023','22-Sep-2023', '26-Sep-2023', '27-Sep-2023', '29-Sep-2023', '01-Oct-2023'...
    %'02-Oct-2023', '04-Oct-2023','05-Oct-2023', '06-Oct-2023', '07-Oct-2023_1', '09-Oct-2023', '10-Oct-2023', '12-Oct-2023','13-Oct-2023'... 
    %'15-Oct-2023', '15-Oct-2023_1','17-Oct-2023','19-Oct-2023', '20-Oct-2023' '23-Oct-2023' '25-Oct-2023' '26-Oct-2023'};

% expID = {'05-Oct-2023', '06-Oct-2023', '07-Oct-2023_1', '09-Oct-2023', '10-Oct-2023', '12-Oct-2023','13-Oct-2023'... 
    % '15-Oct-2023', '15-Oct-2023_1','17-Oct-2023','19-Oct-2023', '20-Oct-2023' '23-Oct-2023' '25-Oct-2023' '26-Oct-2023'};
% 
% expID = {'09-May-2023','11-May-2023','12-May-2023'};
expID = {'14-Aug-2023' '15-Aug-2023' '17-Aug-2023' '21-Aug-2023' '22-Aug-2023' '24-Aug-2023' '28-Aug-2023' '30-Aug-2023' '31-Aug-2023' '01-Sep-2023' '03-Sep-2023' '04-Sep-2023'...
    '14-Sep-2023' '26-Sep-2023' '27-Sep-2023' '03-Oct-2023' '04-Oct-2023' '17-Oct-2023' '19-Oct-2023' '24-Oct-2023'};
lastDay ='24-Oct-2023';

hva ='';
flag = 0 ; %whether or not you need to compile (else - load)
lickTrig = 1;
%% individual animal analysis
if ~lickTrig 
    for i = 1:length(Animals)
        Animal = Animals{i};

        %     prev_path = [fPath filesep Animal filesep 'behavior']; % EK Jan 24
        %     list = dir(prev_path); % changed to dir from ls
        %     explist = []; %list of all dates for the animal (not just new ones given as input)
        %     for i = 1:height(list)
        %         if double(list(i,1).name(1)) >= 48 && double(list(i,1).name(1) <=57) %makes sure it's only the dates
        %             explist = [explist;{list(i,1).name}];
        %         end
        %     end
        %     if ~strcmp(explist{end},lastDay) %cut off list at given last day
        %         explist(end) = [];
        %     end
        % expID = explist';
        % plotDffTracesAndPeaks_beh('hit', Animal, expID(i), hva, 0);
        plotDffTracesAndPeaks_beh('hit', Animal, expID, hva, 0);
        plotDffTracesAndPeaks_beh('miss', Animal, expID, hva, 0);
        plotDffTracesAndPeaks_beh('FA', Animal, expID, hva, 0);
        plotDffTracesAndPeaks_beh('CR', Animal, expID, hva, 0);
        plotDffTracesAndPeaks_beh('LL', Animal, expID, hva, 0);
    end
    %%
else
    for i = 1:length(Animals)
        Animal = Animals{i};
        plotDffTracesAndPeaks_beh_RTV2('hit', Animal, expID, hva, 0);
        % plotDffTracesAndPeaks_beh('miss', Animal, expID, hva, 0);
        plotDffTracesAndPeaks_beh_RTV2('FA', Animal, expID, hva, 0);
        % plotDffTracesAndPeaks_beh('CR', Animal, expID, hva, 0);
        plotDffTracesAndPeaks_beh_RTV2('LL', Animal, expID, hva, 0);
    end
end
%% mean df/f traces across ANIMALS lick triggered
% if lickTrig
% cohort = {'M230220_1'; 'M230209_1'};
% type = 'hit'; typeo = 'FA';
% expIDs = {'20-Oct-2023';'12-May-2023_1'}; %last day for each animal
% hva = [];
% rt = 0;
% 
% temp_gr = cell(14,2); % exp x location
% for a = 1:height(cohort)
%     Animal = cohort{a,:};
%     expID = expIDs{a,:};
%     fPath = 'Y:\haider\Data\analyzedData\EKK\WFI';
%     prev_path = [fPath filesep Animal filesep 'behavior'];
%     % list = cd(prev_path);
%     list = dir(prev_path); % changed to dir from ls
%     gratings = 2; % 0 if barmapping, 1 if passive grating, 2 if active (behavior)
%     nrLocs = 2; % 2 for gratings (binoc/monoc)
% 
%     %load previously compiled data or add on new dates if necessary
%     if rt
%         if isempty(hva)
%             load([prev_path '\' expID '\compiled analyzed data_V1_RT.mat'])
%         else
%             load([prev_path '\' expID '\compiled analyzed data_' hva '_RT.mat'])
%         end
%     else
%         if isempty(hva)
%             load([prev_path '\' expID '\compiled analyzed data_V1.mat'])
%         else
%             load([prev_path '\' expID '\compiled analyzed data_' hva '.mat'])
%         end
%     end
% 
%     comp = analyzed_comp;
%     newExpNr = height(comp);
% 
%     nrSessions = zeros(height(analyzed_comp),1);
%     for exp = 1: height(analyzed_comp)
% 
%         for sessNr = 1:width(analyzed_comp)
%             if ~isempty(analyzed_comp(exp, sessNr).hit)
%                 nrSessions(exp) = nrSessions(exp) + 1;
%             end
%         end
%     end
% 
%     for exp = 1:size(comp,1)
%         for s = 1: size(comp,2)
%             if ~isempty(comp(exp,s).hit)
%                 cont(exp) = length(comp(exp,s).hit.contrasts)/2; % assumes monoc cont nr = binoc cont nr in a session
%                 idx = find(comp(exp,s).hit.contrasts == 0, 1, 'last');
%                 bcontrasts{exp,s} = comp(exp,s).hit.contrasts(1:idx-1);
%                 mcontrasts{exp,s} = comp(exp,s).hit.contrasts(idx:end);
%             end
%         end
%     end
%     b_nrcont = cellfun(@(x) max(length(x)), bcontrasts);
%     m_nrcont = cellfun(@(x) max(length(x)), mcontrasts);
% 
%     data = cell(size(comp,1),nrLocs); % exp x location
%     for exp = 1: size(comp,1) % experiments
%         for s = 1: size(comp,2)% sessions
% 
%             if gratings == 2 % sort binoc and monoc to cell
%                 if ~isempty(comp(exp,s).(type))
%                     if b_nrcont(exp,s) ~= max(max(b_nrcont)) % fill NAN for full contrast case
%                         idxb = unique(cell2mat(bcontrasts(max(max(b_nrcont))== b_nrcont)))';
%                         idxb1 = unique(cell2mat(bcontrasts(min(b_nrcont(b_nrcont>0))== b_nrcont)));
%                         idx = ismember(idxb, idxb1);
%                         comp(exp,s).(type).avgTempTrace.binoc(b_nrcont(exp,s)+1: max(max(b_nrcont)), :) = NaN;
%                         comp(exp,s).(type).avgTempTrace.binoc(idx, :) = comp(exp, s).(type).avgTempTrace.binoc(1:b_nrcont(exp,s),:);
%                         comp(exp,s).(type).avgTempTrace.binoc(idx==0, :) = NaN;
%                     end
%                     if m_nrcont(exp,s) ~= max(max(m_nrcont)) % fill NAN for full contrast case
%                         idxm = unique(cell2mat(mcontrasts(max(max(m_nrcont))== m_nrcont)))';
%                         idxm1 = unique(cell2mat(mcontrasts(min(m_nrcont(m_nrcont>0))== m_nrcont)));
%                         idx = ismember(idxb, idxb1);
%                         comp(exp,s).(type).avgTempTrace.monoc(m_nrcont(exp,s)+1: max(max(m_nrcont)), :) = NaN;
%                         comp(exp,s).(type).avgTempTrace.monoc(idx, :) = comp(exp, s).(type).avgTempTrace.monoc(1:m_nrcont(exp,s),:);
%                         comp(exp,s).(type).avgTempTrace.monoc(idx==0, :) = NaN;
%                     end
%                     %3d cell array; cont x timpoints x session; cell size: exp x location
%                     if ~isempty(comp(exp,s).(type))
%                         if a == 1
%                             aTTb = comp(exp,s).(type).avgTempTrace.binoc;
%                             aTTb(1,:) = comp(exp,s).(typeo).avgTempTrace.binoc(1,:); %add FA
%                             temp_gr{exp,1} = cat(3, temp_gr{exp,1}, aTTb); % concatenate across sessions
% 
%                             aTTm = comp(exp,s).(type).avgTempTrace.monoc;
%                             aTTm(1,:) = comp(exp,s).(typeo).avgTempTrace.monoc(1,:); %add FA
%                             temp_gr{exp,2} = cat(3, temp_gr{exp,2}, aTTm);
%                         else
%                             aTTb = comp(exp,s).(type).avgTempTrace.binoc;
%                             aTTb(1,:) = comp(exp,s).(typeo).avgTempTrace.binoc(1,:); %add FA
%                             temp_gr{exp+11,1} = cat(3, temp_gr{exp+11,1}, aTTb); % concatenate across sessions
% 
%                             aTTm = comp(exp,s).(type).avgTempTrace.monoc;
%                             aTTm(1,:) = comp(exp,s).(typeo).avgTempTrace.monoc(1,:); %add FA
%                             temp_gr{exp+11,2} = cat(3, temp_gr{exp+11,2}, aTTm);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% temp = temp_gr; %both animals! contrasts X timescale X sessions
% 
% tempBD = []; %binoc
% for tt = 1:height(temp)
%     [~,~,h] = size(temp{tt,1});
%     base = temp{tt,1};
%     for ttt = 1:h
%         tempBD = cat(3,tempBD, base(:,:,ttt)); %concatenate days in the third dimension)
%     end
% end
% tempBD = {tempBD};
% 
% tempMD = []; %monoc
% for tt = 1:height(temp)
%     [~,~,h] = size(temp{tt,2});
%     base = temp{tt,2};
%     for ttt = 1:h
%         tempMD = cat(3,tempMD, base(:,:,ttt)); %concatenate days in the third dimension
%     end
% end
% tempMD = {tempMD};
% 
% f1 = figure;
% if isempty(hva)
%     subplot 121
%     [traceb,stdb] = getAvgTraceAcrossSess(temp(exp,1), idxb, 30, 0);
%     title(['Binocular ' type ' - V1']);
%     xlabel('relative time to first lick (s)')
%     ylim([-1.5 4]); hold on;
%     subplot 122
%     [tracem,stdm] = getAvgTraceAcrossSess(temp(exp,2), idxm, 30, 0);
%     title(['Monocular ' type ' - V1']);
%     xlabel('relative time to first lick (s)')
%     ylim([-1.5 4]);
% else
%     subplot 121
%     [traceb,stdb] = getAvgTraceAcrossSess(temp(exp,1), idxb, 30, 0);
%     title([ hva(1:2)]);
%     xlabel('relative time to first lick (s)')
%     ylim([-1.5 4]);
%     subplot 122
%     [tracem,stdm] = getAvgTraceAcrossSess(temp(exp,2), idxm, 30, 0);
%     title([ hva(4:5)]);
%     xlabel('relative time to first lick (s)')
%     ylim([-1.5 4]);
% end
% 
% sgtitle (['Average ' type ' and ' typeo ' Response across animals'])
% cd(['Y:\haider\Data\analyzedData\EKK\WFI\stats\behavior'])
% savefig(f1, ['meandffTrace_animals_' type '_and_' typeo '_RT_' hva '.fig'])
% saveas(f1, ['meandffTrace_animals_' type '_and_' typeo '_RT_' hva '.png'])
% end