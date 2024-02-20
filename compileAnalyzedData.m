function [analyzed_comp, nrSessions] = compileAnalyzedData_NA (fPath, Animal, expID, gratings, hva, flag) 
% compile analyzed mat. file from each session and compile them with different experiments
%INPUTS: (runs from plotActivityVsContrast_beh.m)
% fPath = filepath (i.e. 'Y:\haider\Data\analyzedData\EKK\WFI')
% animal = name of animal to collect data from (i.e. 'M230220_1')
% expID = list of all ADDITIONAL dates to be compiled, cell array (i.e.{'24-May-2023', '25-May-2023'})
% gratings = type of day's (i.e. '2' = behavior)
% EK May23, NA May23
% removed unset third output NA 12/4/23
addpath(fPath)

%load any previously compiled data
if gratings == 1
    prev_path = [fPath filesep Animal filesep 'gratings'];
elseif gratings == 0
    prev_path = [fpath filesep Animal filesep 'barmapping'];
elseif gratings ==2
    prev_path = [fPath filesep Animal filesep 'behavior'];
end

%list = dir(prev_path);
%list = list(3:end);  % remove . and .. % EK

%dateStrings = regexp({list.name}, '(\d{2}-[A-Za-z]{3}-\d{4})', 'match'); % Extract date strings
%[~, sortedIndices] = sort(cellfun(@(x) datenum(x, 'dd-mmm-yyyy'), dateStrings, 'UnifrormOutput', false)); % Sort the list based on datenum values
%list = list(sortedIndices); % Reorder the list based on sortedIndices

%for i = 1: height(list)
 %   if strcmp(list(i).name, expID{1})
  %      startIdx = i;
  %  end
%end

if flag == 1
    if isempty(hva)
        afile = [prev_path filesep expID{end} filesep 'compiled analyzed data.mat']; %EK changed
    else
        afile = [prev_path filesep expID{end} filesep 'compiled analyzed data_' hva '.mat']; %EK changed
    end
else
    if isempty(hva)
        afile = [prev_path filesep expID{end} filesep 'compiled analyzed data.mat']; %EK changed
    else
        afile = [prev_path filesep expID{end} filesep 'compiled analyzed data_' hva '.mat']; %EK changed
    end
end

if exist(afile)
    load(afile)
    analyzed_comp_old = analyzed_comp; %so data isn't written over
    clear analyzed_comp
else % first time compiling data  
    analyzed_comp_old =[];
end

for exp = 1:length(expID)
    if gratings == 1
        dir_path = [fPath filesep Animal filesep 'gratings' filesep expID{exp}];
    elseif gratings == 0
        dir_path = [fPath filesep Animal filesep 'barmapping' filesep expID{exp}];
    elseif gratings == 2
        dir_path = [fPath filesep Animal filesep 'behavior' filesep expID{exp}];
    end

    % get a list of all subfolders
    sessions = dir(dir_path);
    sessions = sessions([sessions(:).isdir]); % keep only directories
    sessions = sessions(3:end); % remove '.' and '..'
    s = [];
    count = 0;
    for i = 1:length(sessions)
        if regexp(sessions(i).name, '^\d+$')
            temp = i; % take the session whose name is an integer only
            sessNr = sessions(temp).name;
            sessNr = str2double(sessNr);
            
            if isempty(hva)
                p = [dir_path filesep num2str(sessNr)];
            else
                p = [dir_path filesep num2str(sessNr) filesep hva];
            end
            if isfolder(p)
            cd(p)
            if gratings == 1 %passive
                cfile = ['analyzed_passive_' num2str(sessNr) '.mat'];
                if exist(cfile)
                    analyzed_comp(exp, sessNr)= load(cfile);
                end
            elseif gratings == 0 %barmapping
                cfile = ['analyzed_barmapping_' num2str(sessNr) '.mat'];
                if exist(cfile)
                    analyzed_comp(exp, sessNr)= load(cfile, 'avgStimResponse','avgTempTrace','cont1','roi');
                end
            elseif gratings ==2 %behavior
                if isempty(hva)
                    cfile = ['analyzed_active_' num2str(sessNr) '.mat'];
                    if isfile(cfile) || isfile([cfile(1:end-4) '_whole.mat'])
                        if isfile(cfile)
                            analyzed_comp(exp, sessNr).whole= load((cfile));
                        else 
                            analyzed_comp(exp, sessNr).whole= load([cfile(1:end-4) '_whole.mat']);
                        end 
                        analyzed_comp(exp, sessNr).hit= load(['analyzed_active_' num2str(sessNr) '_hit.mat']);
                        analyzed_comp(exp, sessNr).miss= load(['analyzed_active_' num2str(sessNr) '_miss.mat']);
                        analyzed_comp(exp, sessNr).CR= load(['analyzed_active_' num2str(sessNr) '_CR.mat']);
                        analyzed_comp(exp, sessNr).LL= load(['analyzed_active_' num2str(sessNr) '_Late.mat']);
                        analyzed_comp(exp, sessNr).FA= load(['analyzed_active_' num2str(sessNr) '_FA.mat']);
                    end
                else
                    cfile = ['analyzed_active_' num2str(sessNr) '_' hva '.mat'];
                     if exist(cfile) || exist([cfile(1:end-4) '_whole.mat']) %edit for hva; EK oct23
                         if isfile(cfile)
                             analyzed_comp(exp, sessNr).whole= load((cfile));
                         else
                             analyzed_comp(exp, sessNr).whole= load([cfile(1:end-4) '_whole.mat']);
                         end
                        analyzed_comp(exp, sessNr).whole= load(cfile);
                        analyzed_comp(exp, sessNr).hit= load(['analyzed_active_' num2str(sessNr) '_' hva '_hit.mat']);
                        analyzed_comp(exp, sessNr).miss= load(['analyzed_active_' num2str(sessNr) '_' hva '_miss.mat']);
                        analyzed_comp(exp, sessNr).CR= load(['analyzed_active_' num2str(sessNr) '_' hva '_CR.mat']);
                        analyzed_comp(exp, sessNr).LL= load(['analyzed_active_' num2str(sessNr) '_' hva '_Late.mat']);
                        analyzed_comp(exp, sessNr).FA= load(['analyzed_active_' num2str(sessNr) '_' hva '_FA.mat']);
                    end
                end
            end
            end 
            count = count +1;
            s = [s sessNr];
        end
        sessNr_comp{exp} = s;
    end

end

if exist(afile)
    %if the number of sessions are different (i.e. width of structures don't match
    if width(analyzed_comp) < width(analyzed_comp_old)
        for exp = 1:length(expID)
            for sessNr = (width(analyzed_comp) + 1):width(analyzed_comp_old)
                analyzed_comp(exp, sessNr).whole= [];
                analyzed_comp(exp, sessNr).hit= [];
                analyzed_comp(exp, sessNr).miss= [];
                analyzed_comp(exp, sessNr).CR= [];
                analyzed_comp(exp, sessNr).LL= [];
                analyzed_comp(exp, sessNr).FA= [];
            end
        end
    elseif width(analyzed_comp) > width(analyzed_comp_old)
        for exp = 1:length(expID)
            for sessNr = (length(analyzed_comp_old) + 1):length(analyzed_comp)
                analyzed_comp_old(exp, sessNr).whole= [];
                analyzed_comp_old(exp, sessNr).hit= [];
                analyzed_comp_old(exp, sessNr).miss= [];
                analyzed_comp_old(exp, sessNr).CR= [];
                analyzed_comp_old(exp, sessNr).LL= [];
                analyzed_comp_old(exp, sessNr).FA= [];
            end
        end
    end
    analyzed_comp = [analyzed_comp_old; analyzed_comp];
end

nrSessions = zeros(height(analyzed_comp),1);
for exp = 1: height(analyzed_comp)
    for sessNr = 1:width(analyzed_comp)
         if ~isempty(analyzed_comp(exp, sessNr).whole)
            nrSessions(exp) = nrSessions(exp) + 1;
        end
    end
end 

if isempty(hva)
    save ([dir_path filesep 'compiled analyzed data.mat'], 'analyzed_comp');
else 
    save ([dir_path filesep 'compiled analyzed data_' hva '.mat'], 'analyzed_comp');
end

