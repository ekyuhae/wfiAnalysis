%% WFI barmapping analysis 
% 1) load the imaging data and preprocess the data. (hemodynamic correction,
% motion correction, spatial downsample)
% 2) calculate dF/F 
% 3) get average activity colormap and plot average dF/F traces 
%    optional: fourier analysis of dF/F traces 
% 4) open analysis GUI and can perform ROI analysis and make movies.

% updated version from basicAnalysis code. 
% EK Feb23 
%% PART I
close all; clear; clc
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

Animal = 'M221118_1';
expID = '02-Feb-2023_1';
cPath = [Animal filesep 'barmapping' filesep expID];
img.savePath = ['Y:\haider\Data\analyzedData\EKK\WFI' filesep cPath]; % directory where you save analysis result
img.fPath = ['Y:\haider\Data\WFI\Mapping\Animals' filesep cPath]; % directory for WFI dataset 
img.fName = 'Frames_2_399_309_uint16'; %name of imaging data files.
img.stimLine = 5; %analog line that contains stimulus trigger.
img.trigLine = [2 3]; %analog lines for blue and violet light triggers.
img.duration = 780;
img.plotChans = true; %flag to show separate channels when loading dual-wavelength data in each trial
img.sRate = 30; % frame rate in Hz
img.downSample = 1; %spatial downsampling factor
img.hemoCorrect = true; %hemodynamic correction is optional (this only works with dual-color data in raw datasets).
img.fileExt = '.dat'; %type of video file. Use '.dat' for binary files (also works for .tif or .mj2 files)
img.preProc = false; %case if data is single channel and can be loaded directly (this is only true for the pre-processed example dataset).
img.photobleaching = false; % do photobleaching assessment (optional)
img.snapshot = true; % true if you want to open snapshot
img.handles = false;  % true if you want to open WFI settings 
plotAnalog = true; 

if img.snapshot
    imshow(imread([img.fPath filesep 'Snapshot_1.jpg']));
end 
if img. handles 
    load([img.fPath filesep 'handles.mat']);
end 
%%
expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal '\2023-02-02\1'];        
sessID = ['2023-02-02_1_' Animal];

load ([expPath filesep sessID '_Timeline.mat']);
load ([expPath filesep sessID '_Block.mat']);
load ([expPath filesep sessID '_parameters.mat']);
%% check photobleaching effect
if img.photobleaching 
    [~, allTimes] = checkSessionBleaching_EK(img.fPath, Animal, img.savePath); 
end 
%% load imaging data
tic 
rawFiles = dir([img.fPath filesep img.fName '*']); %find data files
load([img.fPath filesep 'frameTimes_0001.mat'], 'frameTimes','imgSize') %get size of imaging data
dataSize = floor(imgSize ./ img.downSample); %adjust for downsampling

nrTrials = length(rawFiles); %nr of trials

nrFrames = img.duration * img.sRate; %frames per trial
allData = NaN(dataSize(1),dataSize(2),nrFrames, nrTrials, 'single'); %pre-allocate data array for all trials
for trialNr = 1:nrTrials
    
    aFile = [img.fPath filesep 'Analog_' num2str(trialNr) '.dat']; %current file to be read
    [~,Analog] = loadRawData(aFile,'Analog'); %load analog data
    Analog = double(Analog);
    % baselinePix = [267 198];
    % stimOn = (img.preStim*img.sRate); %frames before stimulus onset
     
    if plotAnalog
    f = figure; plot (Analog(1,:), Analog(2:end, :))
    xlabel('Time (s)'); ylabel('Voltage (mV)')
    legend ('blue LED', 'violet LED', 'cam strobe', 'photodetector')
    end
    [timepoints, stimOn] = getPhotoStartFixed_WFI(img, Analog); %frames before MC starts
    
    [~,~,a] = fileparts(rawFiles(trialNr).name); %get data type (should be .dat or .tif for raw. also works with .mj2 for compressed data)
       
    if ~img.preProc %no preprocessing. assuming two wavelengths for expsure light (blue and violet).       
        
        [bData,bFrameTimes,vData,vFrameTimes] = splitChannels_EK(img,trialNr,a, Analog);  
        [bData, vData] = motionCorrect(bData, vData); %perform motion correction for both channels        
        
        if img.hemoCorrect 
            baseline = floor(length(bData)/2):length(bData);
% %             for i = 1:length(block.trial)
%             %
%             %                 t(i) = block.trial(i).stimulusStartedTime;
%             %             end
%             t1 =[];
%             for i = 1:length(block.trial)
%                 if  block.trial(i).condition.colour(1) ~= 128
%                     t = block.trial(i).stimulusStartedTime;
%                     t1 = [t1;t];
%                 end
%             end
%             if bFrameTimes (1) < vFrameTimes (1) % aligning stimulus start time to frame start time (ms)
%                 time = (t1+stimOn)*1e3 + bFrameTimes(1);
%             else
%                 time =  (t1+stimOn)*1e3 + vFrameTimes(1);
%             end
%             temp = findStimOnFrame(time, 50, bFrameTimes);
%             tmp = [];
%             for i =1:length(temp)
%                 temp1 = temp(i)+1:temp(i)+2;
%                 tmp = [tmp temp1];
%             end
%             idx = 1:length(bData);
%             for i = 1:length(tmp)
%             idx(find(idx == tmp(i))) = [];
%             end 
            data = Widefield_HemoCorrect(bData,vData,baseline, 5); %perform hemodynamic correction for individual pixels
        else % use blue only data
            data = bData;
        end
        
    elseif img.preProc %pre-processed data. simply load all available data and skip motion correction ect.       
        
        cFile = [img.fPath filesep 'frameTimes_' num2str(trialNr, '%04i') '.mat']; %need size of imaging data
        load(cFile, 'imgSize');        
        
        cFile = [img.fPath filesep img.fName '_' num2str(trialNr, '%04i') img.fileExt]; %current file to be read
        [~, data] = loadRawData(cFile, 'Frames', 'uint16', imgSize); %load imaging data         
    else
        error('Could not read number of channels from filename or channelnumber is >2.'); %for this to work filenames should contain the channelnumber after a _ delimiter
    end
    
    %spatially downsample imaging data
    data = arrayResize(data, img.downSample); %this reduces resolution to ~80um / pixel
    if imgSize(end) < nrFrames
        allData(:,:,1:imgSize(end),trialNr) = data;
    else
        allData = data(:,:,1:end);
%         allData(:,:,:,trialNr) = data(:,:,1:nrFrames);
    end
clear data bData vData

if ~img.hemoCorrect %dF/F is automatically computed during hemodynamic correction       
    % compute dF/F by subtracting and dividing the pre-stimulus baseline
    baselineAvg = nanmean(allData(:,:, baselineDur),3);
    allData = reshape(allData, dataSize(1),dataSize(2), []); %merge all frames to subtract and divide baseline
    allData = bsxfun(@minus, allData, baselineAvg); % subtract baseline
    allData = bsxfun(@rdivide, allData, baselineAvg); % divide baseline
    allData = reshape(allData, dataSize(1),dataSize(2),[],nrTrials); %shape back to initial form    
end

sPath = [img.savePath filesep num2str(trialNr)];
mkdir(sPath) 
save([sPath filesep 'preprocessed_WFIdata_' num2str(trialNr,'%04i') '.mat'], 'allData', 'img', 'bFrameTimes', 'vFrameTimes', 'stimOn', 'timepoints')
savefig(f, [sPath filesep 'Analog' num2str(trialNr,'%04i') '.fig'])
fprintf('%d / %d file preprocessed\n', trialNr, nrTrials)
clear allData bFrameTimes vFrameTimes    
end
toc; 
