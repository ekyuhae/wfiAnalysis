%% WFI Behavior Experiment Pipeline 
% during training BEFORE block switch phase which has two locations; phase 0~2 (passive, ecc shift,
% contrast)
% ensure the barmapping analysis has been done beforehand so that you can
% choose ROI based on the barmapping spatial data acquired from the same
% day of behavioral experiment 
% change AnimalID, expID, expID_ret, sessNr, savePath, expInfoNr, and date
% EK Aug23
% for anything before switching

clc; clear; close all

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

% change per animal
Animal = 'M230831_1';
expID_ret = '06-Sep-2023'; % for retinotopy registration and identifying the visual areas 

% change per day
expID = '19-Sep-2023'; % this is WFI experiment ID format

% change per session
sessNr = 1; % enter the imaging session number(X) from Frames_2_xxx_xxx_uint16_000X 
expInfoNr = '1'; % enter the corresponding expInfo session number (does not always match with WFI sessNr) 
stimLocation = 70; %CHANGE

date = datestr(expID(1:11), 'yyyy-mm-dd'); % expInfo date
pupilAnalysis = false; % true for aligning puil trace to image data
stimulusTriggered = true;
reactionTriggered = false;
hva = RL_AM; %higher visual areas being analyzed (i.e. 'RL_AM' or 'LM_PM' or '' for V1)
fwindow = 4; % response window for taking avg spatial map
flag = 1; % 1 if imaged only visual areas
directo ='Y:\haider\Data\analyzedData\EKK\WFI\';
dayType = 'active';

%% create data saving folder and load the wfi/behavior data
cPath = [Animal filesep 'behavior' filesep expID filesep num2str(sessNr)];
if ~isempty(hva)
    if stimulusTriggered
        savePath = ([directo cPath filesep hva]); % enter path to saving the analyzed data
        if ~isfolder(savePath)
            mkdir(savePath)
        end
    end
    if reactionTriggered
        savePathRT = ([directo cPath filesep hva '_RT']);
        if ~isfolder(savePathRT)
            mkdir(savePathRT)
        end
    end
else 
    if stimulusTriggered
        savePath = ([directo cPath]); % enter path to saving the analyzed data
    end
    if reactionTriggered
        savePathRT = ([directo cPath filesep 'V1_RT']);
    end
end 
load([directo cPath filesep 'preprocessed_WFIdata_' num2str(sessNr, '%04i') '.mat'])

% expInfo data
try
catch
    if expID(1:2) ~= date(end-1:end)
        disp('WFI experiment date does not match with expInfo date') % add line to stop running the code
    end
end
expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal filesep date filesep expInfoNr];        
sessID = [date '_' expInfoNr '_' Animal];

%load ([expPath filesep sessID '_Timeline.mat']); % needed for pupil analysis 
load ([expPath filesep sessID '_Block.mat']);
load ([expPath filesep sessID '_parameters.mat']);
[out, beh_all, activeOut] = analyze_daily_JA(Animal, date, dayType, {expInfoNr}); % 8/29 EK changed to attention code from regular analysis code
close all

%% make sure the photodiode onset is correct
% comment this if you preprocess your image data with wfiProcessImage
% instead of barmappingPreproc
% aFile = [img.fPath filesep 'Analog_' num2str(sessNr) '.dat']; %current file to be read
% [~,Analog] = loadRawData(aFile,'Analog'); %load analog data
% Analog = double(Analog);
% 
% [~, stimOn1] = getPhotoStartFixed_WFI(img, Analog); %frames before MC starts
% if stimOn ~= stimOn1
%     a = figure; plot (Analog(1,:), Analog(2:5, :))
%     xlabel('Time (s)'); ylabel('Voltage (mV)')
%     legend ('blue LED', 'violet LED', 'cam strobe', 'photodetector')
%     stimOn_fixed = stimOn1;
% else 
%     stimOn_fixed =stimOn;
% end 
%% align pupil and plot traces together
if pupilAnalysis
     [pupil,me] = pupilAlignWFI(timepoints, stimOn_fixed, Animal, Timeline, sessID(1:10), sessID(12), 0);
   
     allData1 = reshape(allData,[],size(allData,3));
     normpupil = pupil.area/nanmean(pupil.area);
     normME = me.vector/nanmean(me.vector);
     
     s = 0:1/15:size(allData1,2)/15;
     figure;
     subplot(2,1,1)
     plot(s(1:end-1), 100*mean(allData1,1)); hold on;
     t =  0:1/30:size(pupil.time,1)/30;
     plot(t(1:length(normpupil)), normpupil'+10 ); plot(t(1:length(normME)), normME'+20); hold off
     title('before aligned')
     xlabel ('time(s)'); legend('dF/F', 'normalized pupil', 'normalized ME'); 
     legend('boxoff'); legend('orientation','horizontal')
     
     subplot(2,1,2)
     plot(s(1:end-1), 100*mean(allData1,1)); hold on;
     plot(pupil.time, normpupil'+10); plot(pupil.time(1:length(normME)), normME'+20); hold off
     title('after aligned')
     xlabel ('time(s)'); legend('dF/F', 'normalized pupil', 'normalized ME'); 
     legend('boxoff'); legend('orientation','horizontal'); 
     
     savefig([savePath filesep 'pupil trace aligned to WFI.fig'])
     saveas(gcf, [savePath filesep 'pupil trace aligned to WFI.png'])
     
     % align to wfi frametimes 
     if exist('bFrameTimes')
         if bFrameTimes(1)<vFrameTimes(1)
             pupil.time = pupil.time + bFrameTimes(1); 
         else
             pupil.time = pupil.time + vFrameTimes(1);
         end
     end  
end 
clear allData1

%% for whole trials; get spatio temporal activity and save the analyzed data 
% get average response map, average response trace over contrast,
% individual traces as well
tic;
stimOn_fixed = stimOn;

%find frames for hit/misses/FA/CR trials and do the same spatiotemporal analysis

if stimulusTriggered

    [~, paramSorted] = behSortParameter_lrn(beh_all,block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'whole', stimLocation);
    [posntime.H, paramSorted_lrn.H] = behSortParameter_lrn(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'hit', stimLocation);
    [posntime.M, paramSorted_lrn.M] = behSortParameter_lrn(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'miss', stimLocation);
    [posntime.FA, paramSorted_lrn.FA] = behSortParameter_lrn(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'FA', stimLocation); % CR and FA should have frames at 0% contrast
    [posntime.CR, paramSorted_lrn.CR] = behSortParameter_lrn(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'CR', stimLocation);
    [posntime.LL, paramSorted_lrn.LL] = behSortParameter_lrn(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'Late', stimLocation);

    if strcmp(dayType, 'passive') % 2 contrast levels (e.g., 0, 0.85, 0.85, 0.85), 1 location
        % whole,
        % find the row index where non zero value starts (e.g. 2nd row for
        % 85% contrast level) - EK added 
        idx = find(cellfun(@(x) find(x~=stimLocation,1,'first'), paramSorted)==2,1,'first'); %CHANGED
        % tmp =unique(reshape(tmp,[],1));
        % paramSorted{1,idx} =[paramSorted{1,paramSorted(cellfun(@(x) find(x~=0,1,'first'), paramSorted)==2)}(1:2) ; tmp];
        if isequal(paramSorted{1,idx}, paramSorted{1,idx+1})
            for i = idx:length(paramSorted)-1
                paramSorted{1,i+1} = [];
            end
        end
        paramSorted = paramSorted(~cellfun('isempty',paramSorted));
        
        % hits
        idx = cellfun(@(x) find(x~=stimLocation,1,'first'), paramSorted_lrn.H, 'UniformOutput', false);
        idx = find(~cellfun('isempty',idx)==1,1,'first');
        if isequal(paramSorted_lrn.H{1,idx}, paramSorted_lrn.H{1,idx+1})
            for i = idx:length(paramSorted_lrn.H)-1
                paramSorted_lrn.H{1,i+1} = [];
            end
        end
        paramSorted_lrn.H = paramSorted_lrn.H(~cellfun('isempty',paramSorted_lrn.H));
        
        % misses
        idx = cellfun(@(x) find(x~=stimLocation,1,'first'), paramSorted_lrn.M, 'UniformOutput', false);
        idx = find(~cellfun('isempty',idx)==1,1,'first');
        if isequal(paramSorted_lrn.M{1,idx}, paramSorted_lrn.M{1,idx+1})
            for i = idx:length(paramSorted_lrn.M)-1
                paramSorted_lrn.M{1,i+1} = [];
            end
        end
        paramSorted_lrn.M = paramSorted_lrn.M(~cellfun('isempty',paramSorted_lrn.M));

        % false alarms
        idx = cellfun(@(x) find(x~=stimLocation,1,'first'), paramSorted_lrn.FA, 'UniformOutput', false);
        idx = find(~cellfun('isempty',idx)==1,1,'first');
        if isequal(paramSorted_lrn.FA{1,idx}, paramSorted_lrn.FA{1,idx+1})
            for i = idx:length(paramSorted_lrn.FA)-1
                paramSorted_lrn.FA{1,i+1} = [];
            end
        end
        paramSorted_lrn.FA = paramSorted_lrn.FA(~cellfun('isempty',paramSorted_lrn.FA));

        % correct rejects
        idx = cellfun(@(x) find(x~=stimLocation,1,'first'), paramSorted_lrn.CR, 'UniformOutput', false);
        idx = find(~cellfun('isempty',idx)==1,1,'first');
        if isequal(paramSorted_lrn.CR{1,idx}, paramSorted_lrn.CR{1,idx+1})
            for i = idx:length(paramSorted_lrn.CR)-1
                paramSorted_lrn.CR{1,i+1} = [];
            end
        end
        paramSorted_lrn.CR = paramSorted_lrn.CR(~cellfun('isempty',paramSorted_lrn.CR));

        % late licks
        idx = cellfun(@(x) find(x~=0,1,'first'), paramSorted_lrn.LL, 'UniformOutput', false);
        idx = find(~cellfun('isempty',idx)==1,1,'first');
        if isequal(paramSorted_lrn.LL{1,idx}, paramSorted_lrn.LL{1,idx+1})
            for i = idx:length(paramSorted_lrn.LL)-1
                paramSorted_lrn.LL{1,i+1} = [];
            end
        end
        paramSorted_lrn.LL = paramSorted_lrn.LL(~cellfun('isempty',paramSorted_lrn.LL));
    end

    behSpatioTempResp_learning(paramSorted, allData, bFrameTimes, savePath, fwindow, img, 'whole', sessNr, expID_ret, hva, flag);
    behSpatioTempResp_learning(paramSorted_lrn.H, allData, bFrameTimes, savePath, fwindow, img, 'hit', sessNr, expID_ret, hva, 0);
    behSpatioTempResp_learning(paramSorted_lrn.M, allData, bFrameTimes, savePath, fwindow, img, 'miss', sessNr, expID_ret, hva, 0);
    behSpatioTempResp_learning(paramSorted_lrn.FA, allData, bFrameTimes, savePath, fwindow, img, 'FA', sessNr, expID_ret, hva, 0);
    behSpatioTempResp_learning(paramSorted_lrn.CR, allData, bFrameTimes, savePath, fwindow, img, 'CR', sessNr, expID_ret, hva, 0);
    behSpatioTempResp_learning(paramSorted_lrn.LL, allData, bFrameTimes, savePath, fwindow, img, 'Late', sessNr, expID_ret, hva, 0);
    toc;
    % compareMovie(allData)
end


if reactionTriggered
    % Reaction Time frame analysis with different contrast levels
    % align frames to reaction time onset for each hit trial and compare responses
    % to different contrast


    [posntime.H_RT, paramSorted_event.H_RT] = behSortParameter_RT(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'hit');
    [posntime.FA_RT, paramSorted_event.FA_RT] = behSortParameter_RT(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'FA');
    [posntime.LL_RT, paramSorted_event.LL_RT] = behSortParameter_RT(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'Late');
    [~, ~, bFrame.H_RT, mFrame.H_RT] = behSpatioTempResp_event(paramSorted_event.H_RT, allData, bFrameTimes, savePathRT, fwindow, img, 'hit', hva, 1);
    [~, ~, bFrame.FA_RT, mFrame.FA_RT] = behSpatioTempResp_event(paramSorted_event.FA_RT, allData, bFrameTimes, savePathRT, fwindow, img, 'FA', hva,  1);
    [~, ~, bFrame.LL_RT, mFrame.LL_RT] = behSpatioTempResp_event(paramSorted_event.LL_RT, allData, bFrameTimes, savePathRT, fwindow, img, 'Late', hva, 1);
end
