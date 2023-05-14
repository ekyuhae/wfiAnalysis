clc; clear; close all

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

Animal = 'M230220_1';
expID = '07-May-2023';
expID_ret = '04-Apr-2023_1'; % for retinotopy registration get the patch 
sessNr = 1;  
cPath = [Animal filesep 'behavior' filesep expID filesep num2str(sessNr)];
savePath = (['Y:\haider\Data\analyzedData\EKK\WFI\' cPath]); % enter path to saving the analyzed data
load(['Y:\haider\Data\analyzedData\EKK\WFI\' cPath filesep 'preprocessed_WFIdata_' num2str(sessNr, '%04i') '.mat'])

expInfoNr = '4';
date = '2023-05-07';
expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal filesep date filesep expInfoNr];        
sessID = [date '_' expInfoNr '_' Animal];

load ([expPath filesep sessID '_Timeline.mat']);
load ([expPath filesep sessID '_Block.mat']);
load ([expPath filesep sessID '_parameters.mat']);
multiCont = true; % true for multi contrast barmapping analysis 
pupilAnalysis = false; % true for aligning puil trace to image data

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

%% process behavior data

%% for whole trials; get spatio temporal activity and save the analyzed data 
% get average response map, average response trace over contrast,
% individual traces as well
fwindow = 6;
[~, paramSorted] = behSortParameter(block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'whole');
[~, ~, binocOnFrame, monocOnFrame] = behSpatioTempResp(paramSorted, allData, bFrameTimes, savePath, fwindow, img,sessNr, expID_ret);

%% find frames for hit/misses/FA/CR trials and do the same spatiotemporal analysis

[posntime, paramSorted_hit] = behSortParameter(block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed, 'hit');

% do analysis on frames corresponding to hit trials only
[~, ~, bFrameH, mFrameH] = behSpatioTempResp_event(paramSorted, allData, bFrameTimes, savePath, fwindow, img, 'hit');

% compareMovie(allData)
