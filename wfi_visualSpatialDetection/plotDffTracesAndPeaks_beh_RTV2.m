function [temp] = plotDffTracesAndPeaks_beh_RTV2(type, Animal, expID, hva, flag)
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
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK\project'))

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
dateStrings = regexp(explist', '(\d{2}-[A-Za-z]{3}-\d{4})', 'match'); % Extract date strings
[~, sortedIndices] = sort(cellfun(@(x) datenum(x, 'dd-mmm-yyyy'), dateStrings)); % Sort the list based on datenum values
explist = explist(sortedIndices); % Reorder the list based on sortedIndices

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
    [comp, nrSessions] = compileAnalyzedData_NA_RT (fPath, Animal, expID, gratings, hva, flag); % changed from expID to explist EK %removed third output NA 12/4/23
end

for exp = 1:size(comp,1)
    for s = 1: size(comp,2)
        if ~isempty(comp(exp,s).hit)
            cont(exp) = length(comp(exp,s).hit.contrasts)/2; % assumes monoc cont nr = binoc cont nr in a session
            if numel(find(comp(exp,s).hit.contrasts == 0)) > 1
                idx = find(comp(exp,s).hit.contrasts == 0, 1, 'last');
                bcontrasts{exp,s} = comp(exp,s).hit.contrasts(1:idx-1);
                mcontrasts{exp,s} = comp(exp,s).hit.contrasts(idx:end);
            else % only one stim location was used (e.g. early learning phase)
                if comp(exp,s).hit.paramSorted{1,1}(1) < 45 % binoc only
                    bcontrasts{exp,s} = comp(exp,s).hit.contrasts;
                    mcontrasts{exp,s} = [];
                else
                    mcontrasts{exp,s} = comp(exp,s).hit.contrasts;
                    bcontrasts{exp,s} = [];
                end
            end
        end
    end 
end
b_nrcont = cellfun(@(x) max(length(x)), bcontrasts);
m_nrcont = cellfun(@(x) max(length(x)), mcontrasts);

idxb = []; idxm = [];
for exp = 1: size(comp,1)
    for s = 1 :size(comp,2)
        idxb = [idxb bcontrasts{exp,s}];
        idxm = [idxm mcontrasts{exp,s}];
    end
end
idxb = unique(idxb);
idxm =unique(idxm);

temp = cell(1,nrLocs); % exp x location
temp_gr = cell(height(comp),nrLocs); % exp x location
data = cell(size(comp,1),nrLocs); % exp x location
for exp = 1: size(comp,1) % experiments
    for s = 1: size(comp,2)% sessions
        if gratings == 2 % sort binoc and monoc to cell
            if ~isempty(comp(exp,s).(type))
                % binoc
                if isfield(comp(exp,s).(type).avgTempTrace, 'binoc') && b_nrcont(exp,s) ~= max(max(b_nrcont)) % fill NAN for full contrast case
                    % idxb = unique(cell2mat(bcontrasts(max(max(b_nrcont))== b_nrcont)))';
                    idxb1 = unique(cell2mat(bcontrasts(min(b_nrcont(b_nrcont>0))== b_nrcont)));
                    idx = ismember(idxb, idxb1);
                    comp(exp,s).(type).avgTempTrace.binoc(b_nrcont(exp,s)+1: max(max(b_nrcont)), :) = NaN;
                    if length(idxb1) == size(comp(exp, s).(type).avgTempTrace.binoc(1:b_nrcont(exp,s),:),1)
                    comp(exp,s).(type).avgTempTrace.binoc(idx, :) = comp(exp, s).(type).avgTempTrace.binoc(1:b_nrcont(exp,s),:);
                    else 
                         comp(exp,s).(type).avgTempTrace.binoc(idx, :) = comp(exp, s).(type).avgTempTrace.binoc(1:length(idxb1),:);
                    end 
                    comp(exp,s).(type).avgTempTrace.binoc(idx==0, :) = NaN;

                elseif isfield(comp(exp,s).(type).avgTempTrace, 'binoc') && any(cellfun(@(x) ~isequal(x, bcontrasts{1}), bcontrasts(2:end)))
                    % idxb = unique(cell2mat(bcontrasts(max(max(b_nrcont))== b_nrcont)))'; % full contrast range that were used across the whole experiments
                    curridxb = bcontrasts{exp,s};
                    idx = ismember(idxb, curridxb); % indexing current experiment/sesssion's contrast range and allocate the trace data to the index
                    comp(exp,s).(type).avgTempTrace.binoc(b_nrcont(exp,s)+1: length(idxb), :) = NaN;
                    if numel(unique(curridxb)) < numel(curridxb) % if it's the case that has same contrast level; e.g. [0 2 2 3 3%]
                        if ~(strcmp(type, 'miss') || strcmp(type, 'hit') || strcmp(type, 'LL'))
                            % comp(exp,s).(type).avgTempTrace.binoc(1, :) = unique(comp(exp, s).(type).avgTempTrace.binoc(1),:);
                        else % if it's hit or miss; exclude 0% row
                            tmp = [zeros(1, length(comp(exp,s).(type).avgTempTrace.binoc)); unique(comp(exp, s).(type).avgTempTrace.binoc(2:b_nrcont(exp,s),:),"rows", 'stable')]; 
                            comp(exp,s).(type).avgTempTrace.binoc(idx==0, :) = NaN;
                            comp(exp,s).(type).avgTempTrace.binoc(idx, :) = tmp; 
                            tmp = [];
                        end
                    else %if it's the case that has full contrast range w not duplicating elements; e.g. [0 5 10 18 33%]
                        comp(exp,s).(type).avgTempTrace.binoc(idx, :) = (comp(exp, s).(type).avgTempTrace.binoc(1:b_nrcont(exp,s),:));
                    end
                    comp(exp,s).(type).avgTempTrace.binoc(idx==0, :) = NaN;
                end

                %monoc 
                if isfield(comp(exp,s).(type).avgTempTrace, 'monoc') && m_nrcont(exp,s) ~= max(max(m_nrcont)) % fill NAN for full contrast case
                    % idxm = unique(cell2mat(mcontrasts(max(max(m_nrcont))== m_nrcont)))';
                    idxm1 = unique(cell2mat(mcontrasts(min(m_nrcont(m_nrcont>0))== m_nrcont)));
                    idx = ismember(idxm, idxm1);
                    comp(exp,s).(type).avgTempTrace.monoc(m_nrcont(exp,s)+1: max(max(m_nrcont)), :) = NaN;
                    if length(idxm1) == size(comp(exp, s).(type).avgTempTrace.monoc(1:m_nrcont(exp,s),:),1)
                        comp(exp,s).(type).avgTempTrace.monoc(idx, :) = comp(exp, s).(type).avgTempTrace.monoc(1:m_nrcont(exp,s),:);
                    else
                        comp(exp,s).(type).avgTempTrace.monoc(idx, :) = comp(exp, s).(type).avgTempTrace.monoc(1:length(idxm1),:);
                    end
                    comp(exp,s).(type).avgTempTrace.monoc(idx==0, :) = NaN;

                elseif isfield(comp(exp,s).(type).avgTempTrace, 'monoc') && any(cellfun(@(x) ~isequal(x, mcontrasts{1}), mcontrasts(2:end)))
                    % idxm = unique(cell2mat(mcontrasts(max(max(m_nrcont))== m_nrcont)))';
                    curridxm = mcontrasts{exp,s};
                    idx = ismember(idxm, curridxm);
                    comp(exp,s).(type).avgTempTrace.monoc(m_nrcont(exp,s)+1: length(idxm), :) = NaN;
                    if numel(unique(curridxm)) < numel(curridxm)
                        if ~(strcmp(type, 'miss') || strcmp(type, 'hit') || strcmp(type, 'LL')) 
                            % comp(exp,s).(type).avgTempTrace.monoc(idx, :) = unique(comp(exp, s).(type).avgTempTrace.monoc(1:m_nrcont(exp,s),:),"rows", 'stable');
                        else
                            tmp = [zeros(1, length(comp(exp,s).(type).avgTempTrace.monoc)); unique(comp(exp, s).(type).avgTempTrace.monoc(2:m_nrcont(exp,s),:),"rows", 'stable')]; 
                            comp(exp,s).(type).avgTempTrace.monoc(idx==0, :) = NaN;
                            comp(exp,s).(type).avgTempTrace.monoc(idx, :) = tmp;
                            tmp = [];
                        end
                    else
                        comp(exp,s).(type).avgTempTrace.monoc(idx, :) = comp(exp, s).(type).avgTempTrace.monoc(1:m_nrcont(exp,s),:);
                    end
                    comp(exp,s).(type).avgTempTrace.monoc(idx==0, :) = NaN;
                end
            end
            %3d cell array; cont x timpoints x session; cell size: exp x location
            if ~isempty(comp(exp,s).(type))
                if isfield(comp(exp,s).(type).avgTempTrace, 'binoc')
                    temp_gr{exp,1} = cat(3, temp_gr{exp,1}, comp(exp,s).(type).avgTempTrace.binoc); % concatenate across sessions
                end
                if isfield(comp(exp,s).(type).avgTempTrace, 'monoc')
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

if ~exist("idxb") || ~exist("idxm")
    idxb = unique(cell2mat(bcontrasts(max(max(b_nrcont))== b_nrcont)))';
    idxm = unique(cell2mat(mcontrasts(max(max(m_nrcont))== m_nrcont)))';
    % if max(idxb) ~= max(idxm)
    %     idxb(end+1) = max(idxm);
    % end
end
% nrContrasts = [max(max(b_nrcont)) max(max(m_nrcont))];
nrContrasts = [length(idxb) length(idxm)];
for exp = 1: size(avg_sessions,1)
    for j = 1: size(avg_sessions, 2) % locations
        if ~isempty(avg_sessions{exp,j})
            for i = 1:nrContrasts(j)
                [peak_value{exp,j}(i,:), peak_index] = getPeakPostStim(avg_sessions{exp,j}(i,:), 15, 1, 2,[]);
                std_sessions_peak{exp,j}(i,:) = std_sessions{exp,j}(i, peak_index);
            end
        end
    end
end
explist = expID; % EK Oct 23
directory = [fPath filesep Animal filesep 'behavior'];
i = 1;
for exp = 1: length(explist)%newExpNr%:length(avg_sessions) %NA 12/4/23
    if any((any(cellfun(@(x) ~isempty(x), peak_value))))

        f = plotDffVsContrast (idxb, idxm, peak_value(exp,:), std_sessions_peak(exp,:), explist(exp));
        
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

for exp = 1: length(avg_sessions)
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
        if ~isempty(temp{exp,i})
        avg_exp{1,i} = cat(3, avg_exp{1,i}, temp{exp,i}); % concatenate over sessions
        end 
    end
end

avg_days = cellfun(@(x) nanmean(x, 3), avg_exp, 'UniformOutput', false);
std_days =  cellfun(@(x) nanstd(x, 0, 3), avg_exp, 'UniformOutput', false);

for j = 1: size(avg_days, 2) % locations
    for i = 1:nrContrasts(j)
        [peak_value_days{j}(i), peak_index] = getPeakPostStim(avg_days{j}(i,:), 15, 1, 2, []);
        std_days_peak{j}(i) = std_days{j}(i, peak_index);
        % sem_days(i,j) = std_days{j}(i, peak_index)/sqrt(length(std_days{j}(i, peak_index)));
    end
end

plotDffVsContrast (idxb, idxm, peak_value_days, std_days_peak, []);
sgtitle(['Average Response Activity Across Days - ' Animal]);
    yLim = [-1 8]; % Define the y-axis limit
    ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
    set(ax,'ylim',yLim);

cd(directory)
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
% for tt = 1:height(temp)
%     [~,~,h] = size(temp{tt,1});
%     base = temp{tt,1};
%     if ~isempty(temp{tt,1})
%         for ttt = 1:h
%             tempBD = cat(3,tempBD, base(:,:,ttt)); %concatenate days in the third dimension)
%         end
%     end
% end
% tempBD = {tempBD};
tempBD = avg_exp(1,1);

% tempMD = []; %monoc
% for tt = 1:height(temp)
%     [~,~,h] = size(temp{tt,2});
%     base = temp{tt,2};
%     if ~isempty(temp{tt,2})
%         for ttt = 1:h
%             tempMD = cat(3,tempMD, base(:,:,ttt)); %concatenate days in the third dimension
%         end
%     end
% end
% tempMD = {tempMD};
tempMD = avg_exp(1,2);

f1 = figure;
if isempty(hva)
    subplot 121
    [traceb,stdb] = getAvgTraceAcrossSess(tempBD, idxb, 30, 0);
    title(['Binocular ' type ' - V1']);
        xlabel('relative time to first lick (s)')

    subplot 122
    [tracem,stdm] = getAvgTraceAcrossSess(tempMD, idxm, 30, 0);
    title(['Monocular ' type ' - V1']);
        xlabel('relative time to first lick (s)')

else
    subplot 121
    [traceb,stdb] = getAvgTraceAcrossSess(tempBD, idxb, 30, 0);
    title(['Binocular ' type ' - ' hva(1:2)]);
        xlabel('relative time to first lick (s)')

    subplot 122
    [tracem,stdm] = getAvgTraceAcrossSess(tempMD, idxm, 30, 0);
    title(['Monocular ' type ' - ' hva(4:5)]);
        xlabel('relative time to first lick (s)')

end
    yLim = [-2.5 8]; % Define the y-axis limit
    ax = findobj(gcf,'type','axes'); % Find all axes handles in the figure
    set(ax,'ylim',yLim);

sgtitle (['Average ' type ' Response over N = ' num2str(length(explist)) ' days (end at ' explist{end} ')'])
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
    save(['compiled ROI traces_days_' type '.mat'], 'avg_exp', 'avg_days', 'std_days', 'expID', 'idxm', 'idxb', 'avg_sessions', 'std_sessions')
else
    if exist(['compiled ROI traces_days_' type '_' hva '.mat'])
        prev = load(['compiled ROI traces_days_' type '_' hva '.mat'], 'expID');
        expID = [prev.expID expID];
    end
    save(['compiled ROI traces_days_' type '_' hva '.mat'], 'avg_exp', 'avg_days', 'std_days', 'expID', 'idxm', 'idxb', 'avg_sessions', 'std_sessions')
end
