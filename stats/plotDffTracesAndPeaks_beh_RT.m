function plotDffTracesAndPeaks_beh_RT(type, Animal, expID, hva, flag)
% Inputs:
% [1] type = 'miss', 'hit', 'FA', 'CR' or 'LL', character vector
% [2] Animal = 'M230130_2', character vector
% [3] expID = {'26-Feb-2023_3', '27-Feb-2023'}, cell array - only need to put in the last date
% [4] hva = 0 or 1. 0 = v1 analysis, 1 = hva analysis (LM, RL, AM, PM)

%Used in:
% plotActivityVsContrast_behV2.m

% EK 23
% NA 5/30/23 - change to type input
% EK 5/31/23 - integrated getAvgTempTrace function 
% NA 6/1/2023 - applied getAvgTrace across days

%% load analyzed individual session's dataset and compile them
fPath = 'Y:\haider\Data\analyzedData\EKK\WFI';
prev_path = [fPath filesep Animal filesep 'behavior'];
% list = cd(prev_path);
list = dir(prev_path); % changed to dir from ls
gratings = 2; % 0 if barmapping, 1 if passive grating, 2 if active (behavior)
nrLocs = 2; % 2 for gratings (binoc/monoc)
%% compile data

explist = []; %list of all dates for the animal (not just new ones given as input)
for i = 1:height(list)
    if double(list(i,1).name(1)) >= 48 && double(list(i,1).name(1) <=57) %makes sure it's only the dates
        explist = [explist;{list(i,1).name}];
    end
end
if ~strcmp(explist{end},expID{end}) %cut off list at given last day
     explist(end) = [];
end

try %load previously compiled data or add on new dates if necessary
    if isempty(hva)
        load([prev_path '\' expID{end} '\compiled analyzed data_V1_RT.mat'])
    else 
        load([prev_path '\' expID{end} '\compiled analyzed data_' hva '_RT.mat'])
    end 
    comp = analyzed_comp;
    newExpNr = height(comp);

    nrSessions = zeros(height(analyzed_comp),1);
    for exp = 1: height(analyzed_comp)
        
        for sessNr = 1:width(analyzed_comp)
            if ~isempty(analyzed_comp(exp, sessNr).hit)
                nrSessions(exp) = nrSessions(exp) + 1;
            end
        end
    end
catch
    [comp, nrSessions, newExpNr] = compileAnalyzedData_NA_RT (fPath, Animal, explist, gratings, hva, flag); % changed from expID to explist EK
end

for exp = 1:size(comp,1)
    for s = 1: size(comp,2)
        if ~isempty(comp(exp,s).hit)
            cont(exp) = length(comp(exp,s).hit.contrasts)/2; % assumes monoc cont nr = binoc cont nr in a session
            idx = find(comp(exp,s).hit.contrasts == 0, 1, 'last');
            bcontrasts{exp,s} = comp(exp,s).hit.contrasts(1:idx-1);
            mcontrasts{exp,s} = comp(exp,s).hit.contrasts(idx:end);
        end
    end
end
b_nrcont = cellfun(@(x) max(length(x)), bcontrasts);
m_nrcont = cellfun(@(x) max(length(x)), mcontrasts);

%%
temp = cell(1,nrLocs); % exp x location
temp_gr = cell(height(comp),nrLocs); % exp x location
data = cell(size(comp,1),nrLocs); % exp x location
for exp = 1: size(comp,1) % experiments
    %%
    for s = 1: size(comp,2)% sessions

        if gratings == 2 % sort binoc and monoc to cell
            if ~isempty(comp(exp,s).(type))
                if b_nrcont(exp,s) ~= max(max(b_nrcont)) % fill NAN for full contrast case
                    idxb = unique(cell2mat(bcontrasts(max(max(b_nrcont))== b_nrcont)))';
                    idxb1 = unique(cell2mat(bcontrasts(min(b_nrcont(b_nrcont>0))== b_nrcont)));
                    idx = ismember(idxb, idxb1);
                    comp(exp,s).(type).avgTempTrace.binoc(b_nrcont(exp,s)+1: max(max(b_nrcont)), :) = NaN;
                    comp(exp,s).(type).avgTempTrace.binoc(idx, :) = comp(exp, s).(type).avgTempTrace.binoc(1:b_nrcont(exp,s),:);
                    comp(exp,s).(type).avgTempTrace.binoc(idx==0, :) = NaN;
                end
                if m_nrcont(exp,s) ~= max(max(m_nrcont)) % fill NAN for full contrast case
                    idxm = unique(cell2mat(mcontrasts(max(max(m_nrcont))== m_nrcont)))';
                    idxm1 = unique(cell2mat(mcontrasts(min(m_nrcont(m_nrcont>0))== m_nrcont)));
                    idx = ismember(idxb, idxb1);
                    comp(exp,s).(type).avgTempTrace.monoc(m_nrcont(exp,s)+1: max(max(m_nrcont)), :) = NaN;
                    comp(exp,s).(type).avgTempTrace.monoc(idx, :) = comp(exp, s).(type).avgTempTrace.monoc(1:m_nrcont(exp,s),:);
                    comp(exp,s).(type).avgTempTrace.monoc(idx==0, :) = NaN;
                end
                %3d cell array; cont x timpoints x session; cell size: exp x location
                if ~isempty(comp(exp,s).(type))
                    temp_gr{exp,1} = cat(3, temp_gr{exp,1}, comp(exp,s).(type).avgTempTrace.binoc); % concatenate across sessions
                    temp_gr{exp,2} = cat(3, temp_gr{exp,2}, comp(exp,s).(type).avgTempTrace.monoc);
                end
            end
        end
    end
end
temp = temp_gr;

%% plots: mean df/f vs contrast across SESSIONS

avg_sessions = cellfun(@(x) nanmean(x, 3), temp, 'UniformOutput', false);
std_sessions =  cellfun(@(x) nanstd(x, 0, 3), temp, 'UniformOutput', false);

% yLim = [-1.5 4];
% ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
% set(ax,'ylim',yLim);


nrContrasts = [max(max(b_nrcont)) max(max(m_nrcont))];
for exp = 1: size(avg_sessions,1)
    for j = 1: size(avg_sessions, 2) % locations
        if ~isempty(avg_sessions{exp})
            for i = 1:nrContrasts(j)
                [peak_value{exp}(i,j), peak_index] = getPeakPostStim(avg_sessions{exp,j}(i,:), 15, 0.5, 2, []);
                std_sessions_peak{exp}(i,j) = std_sessions{exp,j}(i, peak_index);
            end
        end
    end
end

directory = [fPath filesep Animal filesep 'behavior'];
i = 1;
for exp = 1: newExpNr%:length(avg_sessions)
    if ~isempty(peak_value{exp})
        f = plotDffVsContrast (idxb, idxm, peak_value{exp}*100, std_sessions_peak{exp}*100, explist(exp,:)); 
        cd([directory filesep explist{exp}])
        if isempty(hva)
             if ~isfolder('V1_RT_compiled_sessions')
                mkdir('V1_RT_compiled_sessions')
             end
             cd([directory filesep explist{exp} filesep 'V1_RT_compiled_sessions'])
            savefig(f, ['meandff_over_contrasts_sessions_' type '.fig'])
            saveas(f, ['meandff_over_contrasts_sessions_' type '.png'])
        else
            if ~isfolder(hva)
                mkdir([hva '_RT_compiled_sessions'])
            end
            cd([directory filesep explist{exp} filesep hva '_RT_compiled_sessions'])
            savefig(f, ['meandff_over_contrasts_sessions_' type '_' hva '.fig'])
            saveas(f, ['meandff_over_contrasts_sessions_' type '_' hva '.png'])
            i = i+1;
        end
    end
end
%% plots: mean df/f traces across SESSION

for exp = 1: newExpNr %length(avg_sessions)
    if ~isempty(temp{exp})
        f1 = figure;
        if isempty(hva)
            subplot 121
            [traceb,stdb] = getAvgTraceAcrossSess(temp(exp,1), idxb, 30, 0);
            title(['Binocular ' type ' - V1']);
               xlabel('relative time to first lick (s)')
            subplot 122
            [tracem,stdm] = getAvgTraceAcrossSess(temp(exp,2), idxm, 30, 0);
            title(['Monocular ' type ' - V1']);   xlabel('relative time to first lick (s)')
        else
            subplot 121
            [traceb,stdb] = getAvgTraceAcrossSess(temp(exp,1), idxb, 30, 0);
            title(['Binocular ' type ' - ' hva(1:2)]);   xlabel('relative time to first lick (s)')
            subplot 122
            [tracem,stdm] = getAvgTraceAcrossSess(temp(exp,2), idxm, 30, 0);
            title(['Monocular ' type ' - ' hva(4:5)]);   xlabel('relative time to first lick (s)')
        end
        sgtitle (['Average ' type ' Response over N = ' num2str(nrSessions(exp)) ' sessions (' explist{exp} ')'])
        % cd([directory filesep explist{exp}])
        if isempty(hva)
            cd([directory filesep explist{exp} filesep 'V1_RT_compiled_sessions'])
            savefig(f1, ['meandffTrace_sessions_' type '.fig'])
            saveas(f1, ['meandffTrace_sessions_' type '.png'])
        else
            cd([directory filesep explist{exp} filesep hva '_RT_compiled_sessions'])
            savefig(f1, ['meandffTrace_sessions_' type '_' hva '.fig'])
            saveas(f1, ['meandffTrace_sessions_' type '_' hva '.png'])
        end
    end
end
%% plots: mean df/f vs contrast across DAYS

avg_exp = cell(1,size(avg_sessions,2));
for exp =1: size(avg_sessions,1)
    for i = 1: size(avg_sessions,2) % locations
        avg_exp{1,i} = cat(3, avg_exp{1,i}, temp{exp,i}); % concatenate over sessions
    end
end

avg_days = cellfun(@(x) nanmean(x, 3), avg_exp, 'UniformOutput', false);
std_days =  cellfun(@(x) nanstd(x, 0, 3), avg_exp, 'UniformOutput', false);

for j = 1: size(avg_days, 2) % locations
    for i = 1:nrContrasts(j)
        [peak_value_days(i,j), peak_index] = getPeakPostStim(avg_days{j}(i,:), 15, 0.5, 2, []);
        std_days_peak(i,j) = std_days{j}(i, peak_index);
        % sem_days(i,j) = std_days{j}(i, peak_index)/sqrt(length(std_days{j}(i, peak_index)));
    end
end

plotDffVsContrast (idxb, idxm, peak_value_days, std_days_peak, []);
sgtitle(['Average Response Activity Across Days - ' Animal]);
cd(directory)
if isempty(hva)
    if ~isfolder('V1_RT_compiled_days')
        mkdir('V1_RT_compiled_days')
    end
    cd([directory filesep 'V1_RT_compiled_days'])
    savefig(gcf, ['meandff_over_contrasts_days_' type '.fig'])
    saveas(gcf, ['meandff_over_contrasts_days_' type '.png'])
else
    if ~isfolder([hva '_RT_compiled_days'])
        mkdir([hva '_RT_compiled_days'])
    end
    cd([directory filesep hva '_RT_compiled_days'])
    savefig(gcf, ['meandff_over_contrasts_days_' type '_' hva '.fig'])
    saveas(gcf, ['meandff_over_contrasts_days_' type '_' hva '.png'])
end
%% plots: mean df/f traces across DAYS

tempBD = []; %binoc
for tt = 1:height(temp)
    [~,~,h] = size(temp{tt,1});
    base = temp{tt,1};
    for ttt = 1:h
        tempBD = cat(3,tempBD, base(:,:,ttt)); %concatenate days in the third dimension)
    end
end
tempBD = {tempBD};

tempMD = []; %monoc
for tt = 1:height(temp)
    [~,~,h] = size(temp{tt,2});
    base = temp{tt,2};
    for ttt = 1:h
        tempMD = cat(3,tempMD, base(:,:,ttt)); %concatenate days in the third dimension
    end
end
tempMD = {tempMD};

f1 = figure;
if isempty(hva)
    subplot 121
    [traceb,stdb] = getAvgTraceAcrossSess(temp(exp,1), idxb, 30, 0);
    title(['Binocular ' type ' - V1']);
    xlabel('relative time to first lick (s)')
    subplot 122
    [tracem,stdm] = getAvgTraceAcrossSess(temp(exp,2), idxm, 30, 0);
    title(['Monocular ' type ' - V1']);
       xlabel('relative time to first lick (s)')
else
    subplot 121
    [traceb,stdb] = getAvgTraceAcrossSess(temp(exp,1), idxb, 30, 0);
    title(['Binocular ' type ' - ' hva(1:2)]);
       xlabel('relative time to first lick (s)')
    subplot 122
    [tracem,stdm] = getAvgTraceAcrossSess(temp(exp,2), idxm, 30, 0);
    title(['Monocular ' type ' - ' hva(4:5)]);
       xlabel('relative time to first lick (s)')
end
sgtitle (['Average ' type ' Response over N = ' num2str(height(explist)) ' days (end at ' explist{end} ')'])
cd([directory])
if isempty(hva)
    cd([directory filesep 'V1_RT_compiled_days'])
    savefig(f1, ['meandffTrace_days_' type '.fig'])
    saveas(f1, ['meandffTrace_days_' type '.png'])
else
    cd([directory filesep hva '_RT_compiled_days'])
    savefig(f1, ['meandffTrace_days_' type '_' hva '.fig'])
    saveas(f1, ['meandffTrace_days_' type '_' hva '.png'])
end

%%
% yLim = [min(min(peak_value_days))-max(max(std_days_peak)) max(max(peak_value_days))+max(max(std_days_peak))]; % Define the y-axis limit
% ax = findobj(f,'type','axes'); % Find all axes handles in the figure
% set(ax,'ylim',yLim);
if isempty(hva)
    if exist(['compiled ROI traces_days_' type '.mat'])
        prev = load(['compiled ROI traces_days_' type '.mat'], 'expID');
        expID = [prev.expID expID];
    end
    save(['compiled ROI traces_days_' type '.mat'], 'avg_exp', 'avg_days', 'std_days', 'expID', 'idxm', 'idxb')
else
    if exist(['compiled ROI traces_days_' type '_' hva '.mat'])
        prev = load(['compiled ROI traces_days_' type '_' hva '.mat'], 'expID');
        expID = [prev.expID expID];
    end
    save(['compiled ROI traces_days_' type '_' hva '.mat'], 'avg_exp', 'avg_days', 'std_days', 'expID', 'idxm', 'idxb')
end
