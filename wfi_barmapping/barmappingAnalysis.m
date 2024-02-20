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

Animal = 'M221118_1';
expID = '02-Feb-2023_1'; 
expID_ret = '19-Dec-2022_1'; % for retinotopy image registration, enter the retinotopy experiment date so that it can use the retinotopy data
sessNr =1;  
cPath = [Animal filesep 'barmapping' filesep expID filesep num2str(sessNr)];
savePath = (['Y:\haider\Data\analyzedData\EKK\WFI\' cPath]);
if ~isfolder(savePath)
    mkdir(savePath)
end 
load(['Y:\haider\Data\analyzedData\EKK\WFI\' cPath filesep 'preprocessed_WFIdata_' num2str(sessNr, '%04i') '.mat'])
  
expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal '\2023-02-02\1'];        
sessID = ['2023-02-02_1_' Animal];

load ([expPath filesep sessID '_Timeline.mat']);
load ([expPath filesep sessID '_Block.mat']);
load ([expPath filesep sessID '_parameters.mat']);
multiCont = true; % true for multi contrast barmapping analysis (also works for full contrast)
pupilAnalysis = false; % true for aligning puil trace to image data
postStim = 6; % set number of frames for stimulus triggered average map
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

%% finding time stamps of bar condition 
% sort block.trial condition based on each bar position and get time stamps
% of that. call the frame before and after those stimulus onset and get
% average pixel data of that bar location.

close all
[~, paramSorted] = barmapSortParameter(block, parameters, bFrameTimes, vFrameTimes, stimOn_fixed);
[avgStimResponse, avgTempTrace] = barmapMultiCont(paramSorted, allData, bFrameTimes, savePath, 6, img, multiCont, sessNr, expID_ret); %for 1 chan dataset

%% 
% getAvgTraceAcrossStimLocs(avgTempTrace,cont1,30)
% %% show stimulus triggered activity
% preStim = 1:stimOn*(img.sRate/2);
% 
% % colorRange = 0.03; %range of colorscale for dF/F
% stimResponse = cat(3, wAvgStimResponse, bAvgStimResponse);
% avgMap = nanmean(stimResponse(:,:, :),3); %show average activity after stimulus onset
% preavgMap = nanmean(nanmean(allData(:,:, preStim,:),3),4); %show average activity before stimulus onset
% 
% figure; subplot (1,2,1)
% imagesc(preavgMap, [-0.03 0.03]); %show average activity after stimulus onset
% % imageScale(preavgMap, colorRange); 
% colormap(viridis); 
% colorbar
% title('Pre-stimulus average activity')
% axis image
% 
% subplot(1,2,2)
% imagesc(avgMap, [-0.005 0.005]); 
% % imageScale(avgMap, colorRange); colormap(colormap_blueblackred(256)); 
% colorbar; title('Stimulus-triggered average activity'); axis image
% savefig([img.savePath filesep Animal '_activity map.fig']);
% saveas(gcf, [img.savePath filesep Animal '_activity map.png']);
% 
% %%
% figure; 
% %show an activity trace
% meanTrace = squeeze(nanmean(reshape(allData, [], size(allData,3), size(allData,4)),1)); %average from an interesting pixel (this one is for hindpaw area)
% % meanTrace = squeeze(nanmean(reshape(allData,[], size(allData,3), size(allData,4)),1)); %average activity over all pixels
% % timeTrace = ((1:nrFrames) ./ opts.sRate) - opts.preStim; %time in seconds
% % timeTrace = ((1:nrFrames/2) ./ (img.sRate/2)); %time in seconds
% timeTrace = ((1:nrFrames) ./ img.sRate); %time in seconds
% 
% b = length(timeTrace)/5;
% title('activity trace over one pixel')
% plot(timeTrace(1:length(meanTrace)), meanTrace)
% % plotLine = stdshade(meanTrace', 0.5, 'g', timeTrace); hold on %plot average activity
% % hold on 
% % plot([0 0], plotLine.Parent.YLim, '--k'); %show stimulus onset
% % plot(plotLine.Parent.XLim, [0 0], '--k');
% xlabel('time(s)'); ylabel('dF/F'); 
