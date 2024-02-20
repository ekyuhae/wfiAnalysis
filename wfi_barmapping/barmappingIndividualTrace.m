% EK Feb23 
%% load corresponding session from expInfo and do the individual session analysis  
close all; clear; clc
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

Animal = 'M230131_1';
expID = '01-Mar-2023_1';
expID_ret = '23-Feb-2023'; % for retinotopy registration get the patch 
sessNr = 1;  
cPath = [Animal filesep 'barmapping' filesep expID filesep num2str(sessNr)];
savePath = (['Y:\haider\Data\analyzedData\EKK\WFI\' cPath]);
load(['Y:\haider\Data\analyzedData\EKK\WFI\' cPath filesep 'preprocessed_WFIdata_' num2str(sessNr, '%04i') '.mat'])
  
expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal '\2023-03-01\4'];        
sessID = ['2023-03-01_4_' Animal];

load ([expPath filesep sessID '_Timeline.mat']);
load ([expPath filesep sessID '_Block.mat']);
load ([expPath filesep sessID '_parameters.mat']);
multiCont = true; % true for multi contrast barmapping analysis 
pupilAnalysis = false; % true for aligning puil trace to image data

%% make sure the photodiode onset is correct
aFile = [img.fPath filesep 'Analog_' num2str(sessNr) '.dat']; %current file to be read
[~,Analog] = loadRawData(aFile,'Analog'); %load analog data
Analog = double(Analog);
    
[~, stimOn1] = getPhotoStartFixed_WFI(img, Analog); %frames before MC starts
if stimOn ~= stimOn1
    a = figure; plot (Analog(1,:), Analog(2:5, :))
    xlabel('Time (s)'); ylabel('Voltage (mV)')
    legend ('blue LED', 'violet LED', 'cam strobe', 'photodetector')
    stimOn_fixed = stimOn1;
else 
    stimOn_fixed =stimOn;
end 

%%
[~, paramSorted] = barmapSortParameter(block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed);
color = unique(paramSorted(1,:));
for i = 1: floor(length(color)/2)+1
    b = paramSorted(3:end, paramSorted(1,:)==color(i));
    w = paramSorted(3:end, paramSorted(1,:) == color(end-(i-1)));
    if ~isequal(b,w)
        contrasts{i} = vertcat(b,w);
    else 
        contrasts{i} = b;
    end
end 
%% 

cont1 =abs(100*(color - color(color==128))/abs(color(color==128))); % convert color into relative contrast change
cont1 = cont1(1:find(color==128));
% loc = paramSorted(2,find(paramSorted(1,:)==0));
loc = paramSorted(2,1:17); % change accordingly if the number of location is not 17!!!
loc_gr = paramSorted(2,find(paramSorted(1,:)==128));
window = img.sRate*2+1; % if padded as -2s~2s window including stim onset frame
sPath = savePath;

for cont = 1:length(cont1) % get response maps from different contrasts (contrast sorted in descending order)
    if cont1(cont) ~= 0
        for locc = 1:length(loc) % locations
            
            frameTemp =[];
            barFrames = findStimOnFrame(contrasts{cont}(:,locc),50,bFrameTimes);
            for i = 1:length(barFrames) % should be same as number of bar reapeats
                if barFrames(i) ~= 0
                    frameTemp= cat(2,frameTemp, barFrames(i)-img.sRate:barFrames(i)+img.sRate); % for the temporal analysis, save frames for +/- window from stimulus onset
                end
            end
            
            cfile =[sPath(1:end-1) num2str(1) filesep 'analyzed_barmapping_1.mat'];
            % cfile = 'Y:\haider\Data\analyzedData\EKK\WFI\M230130_2\barmapping\09-Mar-2023\1\3ROIs\analyzed_barmapping_s.mat';

            if exist(cfile)
                load(cfile,'roi', 'cord')
            end
            
            temp = reshape(frameTemp', window, [])'; %matrix to get n (number of repeats) traces
            
            idx = find(temp > length(allData));
            if ~isempty(idx)
                temp(end,:) = [];
            end
            
            %             for i = 1:size(temp,2) % avg across number of repeats
            %                 avgTemp(:,:,i) = mean(allData(:,:,temp(:,i)),3);
            %             end
            
            for i = 1:size(temp,1)
                if isempty (roi{locc}) % get mean trace
                    [dFFtrace{cont}(i,:), roi{locc}] = getAvgPixValROI(allData(:,:,temp(i,:)),window, [], cord, sPath); % for first contrast, set roi
                else
                    [dFFtrace{cont,locc}(i,:),~] = getAvgPixValROI (allData(:,:,temp(i,:)), window, roi{locc}, cord, sPath);
                end
            end
            
        end
    else % for zero contrast (control trial)
        
        frameTemp =[];
        barFrames = findStimOnFrame(contrasts{cont},50,bFrameTimes);
        for i = 1:length(barFrames) % should be same as number of bar reapeats
            frameTemp= cat(2,frameTemp, barFrames(i)-img.sRate:barFrames(i)+img.sRate); % for the temporal analysis, save frames for +/- window from stimulus onset
        end
        temp = reshape(frameTemp', window, [])'; %matrix to get n (number of repeats) traces
        
        idx = find(temp > length(allData));
        if ~isempty(idx)
            temp(end,:) = [];
        end
        for i = 1:size(temp,1)
            [dFFtrace{cont,1}(i,:),~] = getAvgPixValROI (allData(:,:,temp(i,:)), window, roi{1}, cord, sPath);
        end
    end
end

save([sPath filesep 'indivTempTraces.mat'], 'dFFtrace')

