%% WFI barmapping analysis 
% 1) load the imaging data and preprocess the data. (hemodynamic correction,
% motion correction, spatial downsample)
% 2) calculate dF/F 
% 3) get average activity colormap and plot average dF/F traces 
%    optional: fourier analysis of dF/F traces 
% 4) open analysis GUI and can perform ROI analysis and make movies.

% updated version from basicAnalysis code. 
% EK Feb23 
%% load corresponding session from expInfo and do the individual session analysis  
close all; clear; clc
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Arvind_Ramesh\AR_function'))


Animal = 'M230131_1';
expID = '02-Mar-2023_1';
expID_ret = '23-Feb-2023'; % for retinotopy registration get the patch 
sessNr = 1;
cPath = [Animal filesep 'gratings' filesep expID filesep num2str(sessNr)];
savePath = (['Y:\haider\Data\analyzedData\EKK\WFI\' cPath]);
load(['Y:\haider\Data\analyzedData\EKK\WFI\' cPath filesep 'preprocessed_WFIdata_' num2str(sessNr, '%04i') '.mat'])

expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal '\2023-03-02\4'];        
sessID = ['2023-03-02_4_' Animal];
load ([expPath filesep sessID '_Timeline.mat']);
load ([expPath filesep sessID '_Block.mat']);
load ([expPath filesep sessID '_parameters.mat']);
multiContrast = true; % true for multi contrast barmapping analysis 
pupilAnalysis = false; % true for aligning puil trace to image data
temporalAnalysis = true;

%% make sure the photodiode onset is correct
aFile = [img.fPath filesep 'Analog_' num2str(sessNr) '.dat']; %current file to be read
[~,Analog] = loadRawData(aFile,'Analog'); %load analog data
Analog = double(Analog);
    
[~, stimOn1] = getPhotoStartFixed_WFI(img, Analog); %frames before MC starts
if stimOn ~= stimOn1
    a = figure; plot (Analog(1,:), Analog(2:end, :))
    xlabel('Time (s)'); ylabel('Voltage (mV)')
    legend ('blue LED', 'violet LED', 'cam strobe', 'photodetector')
    stimOn_fixed = stimOn1;
else 
    stimOn_fixed = stimOn;
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

%% analyze stimulus triggered activity in spatial and temporal space
close all

[~, paramSorted] = barmapSortParameter_gratings(block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed);
[avgStimResponse, avgTempTrace, binocOnFrame, monocOnFrame] = barmapMultiCont_gratings(paramSorted, allData, bFrameTimes, savePath, 6, img,sessNr, expID_ret); % you can also apply full contrast grating
