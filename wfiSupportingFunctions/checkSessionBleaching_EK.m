function [allData, allTimes] = checkSessionBleaching_EK(cPath,Animal,sPath)
%CHECKSESSIONBLEACHING_EK plot mean frame intensity changes throughout session
%cPath {char} absolute path to session folder to source .dat video files
%Animal {char} animal ID (REDUNDANT IN CPATH) used in saved figure filename.
%sPath {char} absolute path to session folder to save figures
%docstring updated Jason Sebek 2023

%% Set basic variables
fclose('all');
fName = 'Frames'; %name format for imaging files
numChans = 2; %number of excitation wavelengths (1 for blue only, 2 for alternating illumination)
dataType = 'uint16'; %type of imaging data usually uint16 or uint8

%% load data
rawVids = dir([cPath filesep fName '_*']); %video files

% get trial numbers for each file and sort in ascending order
trials = zeros(1,length(rawVids));
for x = 1 : length(rawVids)
    [~,a] = fileparts(rawVids(x).name);
    a = textscan(a,'%s','delimiter','_');
    trials(x) = str2double(a{1}{end});
end
[trials,sortIdx] = sort(trials,'ascend');

% get single frame data
allData = cell(1,length(trials));
allTimes = cell(1,length(trials));
for iTrials = 1 : length(trials)
    
    if numChans == 2
        cFile = [cPath filesep rawVids(sortIdx(iTrials)).name];
        try
            load([cPath filesep 'frameTimes_' num2str(trials(iTrials), '%04i')], 'imgSize', 'frameTimes') %get size of imaging file EK 09/22/22 changed 
            [~, cData] = loadRawData(cFile, 'Frames', [], imgSize);
        catch
            [header, cData] = loadRawData(cFile, 'Frames', dataType); %if file was written was old version, the image size is stored in the binary file
            if isempty(header) %if no header is found frametimes should be
                load([cPath filesep 'frameTimes_' num2str(trials(iTrials))], 'frameTimes');
            else
                frameTimes = header(1:end-4);
            end
        end
        if length(size(cData)) == 4
            if size(cData,3) == 3 % convert rgb image to gray
                cData =  squeeze(sum(cat(3, 0.2989 .* cData(:,:,1,:), 0.5870 .* cData(:,:,2,:), 0.1140 .* cData(:,:,3,:)), 3));
            else
                cData = squeeze(cData(:,:,1,:));
            end
        end
        % dont use first and last frames
        cData = cData(:,:,31:end-30);
        frameTimes = frameTimes(31:end-30);
    end
    allData{iTrials} = mean(reshape(cData,[], size(cData,3),1))'; %keep frame average
    allTimes{iTrials} = frameTimes;
end

% combine all trials
allData = cat(1, allData{:});
allTimes = cat(1, allTimes{:});
allTimes = (allTimes - allTimes(1)) .* 1440; %convert to minutes
allData(zscore(diff(allTimes)) > 1) = NaN; %to avoid lines between trials
allData(~isnan(allData)) = smooth(allData(~isnan(allData)), 50); % do some smooth

%% show plot
figure('name', 'check photobleaching')
plot(allTimes,allData, 'linewidth', 2, 'Color', 'k'); 
axis square; xlabel('Time (minutes)'); ylabel('Mean fluorescence');
title(sprintf('Mean frame intensity\ndata path: %s', cPath))
if ~isfolder(sPath)
    mkdir(sPath)
    savefig([sPath filesep Animal 'PhotobleachingAssessment.fig']); %EK
    saveas(gcf, [sPath filesep Animal 'PhotobleachingAssessment.jpg']); %EK
else 
    savefig([sPath filesep Animal 'PhotobleachingAssessment.fig']); %EK
    saveas(gcf, [sPath filesep Animal 'PhotobleachingAssessment.jpg']); %EK
end
