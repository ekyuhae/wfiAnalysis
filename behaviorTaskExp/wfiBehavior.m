%% WFI Behavior Experiment Pipeline 
% ensure the barmapping analysis has been done beforehand so that you can
% choose ROI based on the barmapping spatial data acquired from the same
% day of behavioral experiment 
% EK May23

clc; clear; close all

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

Animal = 'M230220_1';
expID = '15-May-2023'; % this is WFI experiment ID format
expID_ret = '04-Apr-2023_1'; % for retinotopy registration and identifying the visual areas 
% expID_ret = '16-Mar-2023';
sessNr = 5; % enter the imaging session number(X) from Frames_2_xxx_xxx_uint16_000X 
cPath = [Animal filesep 'behavior' filesep expID filesep num2str(sessNr)];
savePath = (['Y:\haider\Data\analyzedData\EKK\WFI\' cPath]); % enter path to saving the analyzed data
load(['Y:\haider\Data\analyzedData\EKK\WFI\' cPath filesep 'preprocessed_WFIdata_' num2str(sessNr, '%04i') '.mat'])

% expInfo data
expInfoNr = '7'; % enter the corresponding expInfo session number (does not always match with WFI sessNr) 
date = '2023-05-15'; 
expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal filesep date filesep expInfoNr];        
sessID = [date '_' expInfoNr '_' Animal];

load ([expPath filesep sessID '_Timeline.mat']);
load ([expPath filesep sessID '_Block.mat']);
load ([expPath filesep sessID '_parameters.mat']);
pupilAnalysis = false; % true for aligning puil trace to image data
fwindow = 6; % response window for taking avg spatial map
[~, totalTrials, ~, beh_all, activeOut] = attention_effects_EK(Animal, date, {expInfoNr});
close all

%% make sure the photodiode onset is correct
% comment this if you preprocess your image data with wfiProcessImage
% instead of barmappingPreproc
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

[~, paramSorted] = behSortParameter([],block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'whole');
[~, ~, binocOnFrame, monocOnFrame] = behSpatioTempResp(paramSorted, allData, bFrameTimes, savePath, fwindow, img,sessNr, expID_ret);

%% find frames for hit/misses/FA/CR trials and do the same spatiotemporal analysis


[posntime.H, paramSorted_event.H] = behSortParameter(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'hit');
[posntime.M, paramSorted_event.M] = behSortParameter(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'miss');
[posntime.FA, paramSorted_event.FA] = behSortParameter(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'FA'); % CR and FA should have frames at 0% contrast
[posntime.CR, paramSorted_event.CR] = behSortParameter(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'CR');
[posntime.LL, paramSorted_event.LL] = behSortParameter(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'Late');

[~, ~, bFrame.H, mFrame.H] = behSpatioTempResp_event(paramSorted_event.H, allData, bFrameTimes, savePath, fwindow, img, 'hit');
[~, ~, bFrame.M, mFrame.M] = behSpatioTempResp_event(paramSorted_event.M, allData, bFrameTimes, savePath, fwindow, img, 'miss');
[~, ~, bFrame.FA, mFrame.FA] = behSpatioTempResp_event(paramSorted_event.FA, allData, bFrameTimes, savePath, fwindow, img, 'FA');
[~, ~, bFrame.CR, mFrame.CR] = behSpatioTempResp_event(paramSorted_event.CR, allData, bFrameTimes, savePath, fwindow, img, 'CR');
[~, ~, bFrame.LL, mFrame.LL] = behSpatioTempResp_event(paramSorted_event.LL, allData, bFrameTimes, savePath, fwindow, img, 'Late');

% compareMovie(allData)
%% Reaction Time frame analysis with different contrast levels
% align frames to reaction time onset for each hit trial and compare responses
% to different contrast

[posntime.H, paramSorted_event.H] = behSortParameter(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'hit');