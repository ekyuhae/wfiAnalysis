%% PhaseMapAnalysis code with integrated WFI system 
% only use it for 2 channel imaging data (modified preprocessing pipeline
% for 2 channel recording version)
% modified spltchannel code to "ret_splitChannels_EK"
% 1) set animal ID, file path (savePath) you want to save your data and the path where
% the original dataset exist (fPath). and change the file name
clc; clear; close all; 
Animal = 'M230220_1';
date = '04-Apr-2023_1';

img.savePath =['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\retinotopy\' date];
if ~isfolder(img.savePath)
    mkdir(img.savePath);
end 
img.fPath = ['Y:\haider\Data\WFI\Mapping\Animals\' Animal '\retinotopy\' date];
img.fName = 'Frames_2_417_296_uint16';
img.plotChans = false;
img.trigLine = [2 3];
img.blueHigh = true;
img.stimLine = 5; %trigger from stimulator
img.alignRes = 10;
img.singleChannel = false;
img.hemoCorrect = true;
img.motionCorr = false;
img.frameRate = 15; %framerate in Hz (this framerate is only for the blue channel!)
img.fileExt = '.dat'; %type of video file. Use '.dat' for binary files (also works for .tif or .mj2 files)
img.checkSessionBleaching = true;
img.makemovie = true;
img.binning = false;
img.refAlign = true;
screenSize = [144.15 80.85]; %size of screen in visual degrees
numCycles = 1;
phaseMapSmth = 2; %smoothing of phasemaps. can help to get better sign maps.
nTrials = 30; %number of trials for each session    
winSize = 5;  %size of the field in mm. This is to determine the spatial binning to get as close to 40pix/mm as possible. This is advised for the visual segmentation code.
rotateImage = 0;  
baselinecorr = true;
%% data source
disp(['Current path: ' img.fPath]);
rawFiles = dir([img.fPath filesep img.fName '*']); %find data files
nrSessions = length(rawFiles); %nr of sessions
trialFrame = cell(nrSessions, nTrials); BaselineFrame = cell(nrSessions, nTrials);
iBarConds1 = [];  
Data = cell(1,nrSessions);
offset = 40; % 40ms for 15Hz fps
%% check for photobleaching 
if img.checkSessionBleaching 
    [data, allTimes] = checkSessionBleaching_EK(img.fPath, Animal, img.savePath); 
end 
    
%% load all dataset and preprocess the imaging data
tic;
for s = 1: nrSessions 
    %% get different bar direction and orentiations for individual trials
    sFile = ls([img.fPath filesep Animal '*settings.mat']);
    StimDataset(1,s) = load([img.fPath filesep sFile(s,:)], 'StimData');
   
    StimDur = StimDataset(s).StimData.VarVals(strcmpi(StimDataset(1).StimData.VarNames,'StimDuration') | strcmpi(StimDataset(1).StimData.VarNames,'trialDuration'),:); %Duration of a given trial
    barFreq = StimDataset(s).StimData.VarVals(strcmpi(StimDataset(1).StimData.VarNames,'cyclesPerSecond'),:); %bar speed in a given trial
    BarOrient = StimDataset(s).StimData.VarVals(strcmpi(StimDataset(s).StimData.VarNames,'BarOrient'),:); %get bar orientation for all trials
    BarDirection = StimDataset(s).StimData.VarVals(strcmpi(StimDataset(s).StimData.VarNames,'BarDirection'),:); %get bar direction for all trials
    
    iBarConds(1,:,s) = BarOrient == 1  & BarDirection == 0; % horizontal moving down
    iBarConds(2,:,s) = BarOrient == 1  & BarDirection == 1; % horizontal moving up
    iBarConds(3,:,s) = BarOrient == 0  & BarDirection == 0; % vertial, moving left to right
    iBarConds(4,:,s) = BarOrient == 0  & BarDirection == 1; % vertical, moving right to left
    iBarConds1 = cat(2, iBarConds1, iBarConds(:,:,s));
    
    Trials = 1 : str2double(StimDataset(1).StimData.handles.NrTrials);
    sRate = round(1 / StimDataset(1).StimData.sRate) / 2; %frame rate in Hz

    %% load dataset
    aFile = [img.fPath filesep 'Analog_' num2str(s) img.fileExt]; %current analog file to be read
    [~,Analog] = loadRawData(aFile,'Analog'); %load analog data
    Session(s).Analog = double(Analog);
    cFile = [img.fPath filesep img.fName '_' num2str(s, '%04i') img.fileExt]; %current imaging file to be read
    
    Session(s).info = load ([img.fPath filesep 'frameTimes_' num2str(s, '%04i') '.mat'], 'imgSize', 'frameTimes'); % get size of imaging data and frame times
    cFile = [img.fPath filesep img.fName '_' num2str(s, '%04i') img.fileExt]; %current file to be read
    [~, Data{1,s}] = loadRawData(cFile, 'Frames', 'uint16', Session(s).info.imgSize); %load imaging data

    frameTimes = Session(s).info.frameTimes*86400*1e3;
    [Session(s).PDTime, stimOn, stimOff] = findStimOns(img, Session(s).Analog);
    if length(stimOff) ~= nTrials % if the frame acquisition stopped earlier than the stimulus offset
        stimOff(nTrials) = Session(s).PDTime(end);
    end 
    Session(s).stimOnTimes = frameTimes (1) + stimOn*1e3; % aligned stimulus timestamps to frametimes timescale 
%         disp (diff(Session(s).stimOnTimes))
    Session(s).stimOffTimes = frameTimes (1) + stimOff*1e3; % aligned stimulus timestamps to frametimes timescale 
    Session(s).trialOnTimes = [frameTimes(1) frameTimes(1) + (stimOff(1:end-1) + str2double(StimDataset(1,s).StimData.handles.ITI))*1e3]; % aligning trials after 1st trial by adding stimoff times + inter trial interval, returns baseline+stim duration
    Session(s).trialOffTimes = Session(s).stimOffTimes;  

    %% identifying trial start frame time (on) and trial end frame time (stim off time)
    for blockTrial = Trials

        [trialFrame{s,blockTrial}, BaselineFrame{s,blockTrial}] = cutFrames_EK(frameTimes, Data, Session, offset, s, blockTrial);

    end 
      %% trial by trial correction
    for blockTrial = Trials       
        if ~ img.singleChannel        
            blueData = cell(1, nTrials); hemoData = cell(1, nTrials);
            [blueData{1,blockTrial},~,hemoData{1,blockTrial},~,~] = ret_splitChannels_EK(img, trialFrame{s,blockTrial}, frameTimes, blockTrial, Session(s).Analog);
            baseline = 1:floor(size(BaselineFrame{s,blockTrial},3)/2);
            if img.hemoCorrect   % perform trial by trial hemodynamic correction %modified by EK  
               
               trialFrame{s,blockTrial} = Widefield_HemoCorrect(blueData{1,blockTrial},hemoData{1,blockTrial},baseline,3);
               clear hemoData blueData
            else
                trialFrame{s,blockTrial} = blueData{1,blockTrial}; clear blueData
                if baselinecorr
                    dataAvg = mean(trialFrame{s,blockTrial}(:,:,baseline),3); %baseline correction
                    trialFrame{s,blockTrial} = bsxfun(@minus, double(trialFrame{s,blockTrial}), dataAvg); % subtract baseline mean
                    trialFrame{s,blockTrial} = bsxfun(@rdivide, double(trialFrame{s,blockTrial}), dataAvg);
                end
            end
        
            rlsize = size(trialFrame{s,blockTrial},3)- length(baseline);
            estsize = StimDur(1)*img.frameRate;
            if rlsize >= estsize
                trialFrame{s,blockTrial}(:,:,baseline) = []; %throw away baseline
                trialFrame{s,blockTrial} = trialFrame{s,blockTrial}(:,:,1:StimDur(1)*img.frameRate);
            else
                diff = estsize-rlsize;
                trialFrame{s,blockTrial}(:,:,baseline(1:end-diff)) = []; %throw away baseline
                trialFrame{s,blockTrial} = trialFrame{s,blockTrial}(:,:,1:StimDur(1)*img.frameRate);
            end
        end

        if img.binning == 1  % perform spatial downsampling
            binSize = floor((max(size(blueData(:,:,1)))/winSize)/40); %compute binsize to get closest to 40 pixels/mm (recommended for segmentation code).
        else
            binSize =1;
        end
%         if winSize > 1 && winSize < inf
%             blueData = arrayResize(blueData,binSize); %do spatial binning
%             hemoData = arrayResize(hemoData,binSize); %do spatial binning
%         end
% 
%         if img.motionCorr   % perform motion correction
%             if s == 1
%                 blueData = single(squeeze(blueData));
%                 blueRef = fft2(median(blueData,3)); %blue reference for motion correction
%                 hemoData = single(squeeze(hemoData));
%                 violetRef = fft2(median(hemoData,3)); %violet reference for motion correction
%             end
%             for iFrames = 1:size(blueData,3)
%                 [~, temp] = dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10); % perform motion correction
%                 blueData(:, :, iFrames) = abs(ifft2(temp));
%                 [~, temp] = dftregistration(violetRef, fft2(hemoData(:, :, iFrames)), 10);
%                 hemoData(:, :, iFrames) = abs(ifft2(temp));
%             end
%             Data{1,s} = blueData;
%         end
    end 
   disp(['Done loading session ' int2str(s) '/' int2str(nrSessions)]); 
end 
toc;

%% 
if img.makemovie % to see simulus evoked response for each bar conditions (optional)
    azimuthAverages = trialAveragedResponse(trialFrame); %outputs numberOfAzimuths x 1 column cell Array
    for i = 1:length(azimuthAverages) % by sweepDirection
        compareMovie(azimuthAverages{i,1}) % make trial-averaged movie for each sweepDirection
    end
end 

disp(['Trials per condition: ' num2str(nTrials) ' - Computing phasemaps every [' num2str(nTrials) '] trials']);

TrialCnts = ones(nrSessions,length(nTrials));
allData = reshape (trialFrame',[],1); 
condCnt = zeros(1,nrSessions); %counter for how many sequences were saved in each condition
iBarConds = logical(iBarConds1); 
clear iBarConds1 s blockTrial stimOnFrame trialOnFrame trialOffFrame frameIdx

%% compute the amount of required frames and collect from data
% for iTrials = 1:length(allData) 
%     if numCycles ~= floor(((1/img.frameRate)*size(allData{iTrials},3)) / (1/barFreq(1))) %make sure that number of cycles matches the current data
%         error('Number of presented cycles does not match the length of available data')
%     end
% end   

%% get average data for each bar conditio and compute fourier transform
a = cell(length(nTrials),1); avgData = cell(nrSessions,1); fTransform = cell(length(nTrials), 1);
for s =1:nrSessions
     trialCnt = 1; 
     for x = Trials             
        a{iBarConds(:,nTrials*(s-1) + trialCnt),1} = allData{nTrials*(s-1) + trialCnt};
        avgData{iBarConds(:,nTrials*(s-1) + trialCnt),1} = cat (4, avgData{iBarConds(:,nTrials*(s-1) + trialCnt),1}, a{iBarConds(:,nTrials*(s-1) + trialCnt),1}); 
        if trialCnt+1 > length(Trials) 
            avgData{iBarConds(:,nTrials*(s-1) + trialCnt),1} = mean(avgData{iBarConds(:,nTrials*(s-1) + trialCnt),1}, 4);
            temp = fft(avgData{iBarConds(:,nTrials*(s-1) + trialCnt),1},[],3);
            fTransform{iBarConds(:,nTrials*(s-1) + trialCnt),1}(TrialCnts(iBarConds(:,nTrials*(s-1) + trialCnt),1),:,:) = squeeze(temp(:,:,numCycles + 1)); 
            clear temp
            TrialCnts(iBarConds(:,nTrials*(s-1) + trialCnt),1) = TrialCnts(iBarConds(:,nTrials*(s-1) + trialCnt),1)+1;
        end 
        trialCnt = trialCnt+1;            
     end       
end

save([img.savePath 'fTransform.mat'],'fTransform'); %keep this to reconstruct phasemaps later if needed
% clear Data

%% do fft analysis to get phase and magnitude maps
for iTrials = 1
    Cnt = 1;
    for iConds = [1 3] %this expects 4 directions to construct horizontal and vertical map (horizontal first)
        for iRuns = 1:size(fTransform{iConds,iTrials},1)
        
           magMaps{Cnt,iRuns} = imrotate(squeeze(abs(fTransform{iConds,iTrials}(iRuns,:,:).*fTransform{iConds+1,iTrials}(iRuns,:,:))),rotateImage); %combined magnitude map.
           if iConds == 3 % azimuth
               a1 = mod(-angle(fTransform{iConds,iTrials}(iRuns,:,:)), 2 * pi); % angle changed from negative to positive get phase angle btw [-pi pi]
               a2 = mod(-angle(fTransform{iConds+1,iTrials}(iRuns,:,:)), 2 * pi);
               a1(isnan(a1)) = 0; a2(isnan(a2)) = 0;
               phaseMaps{Cnt,iRuns} = imrotate(squeeze((a1 - a2) / 2),rotateImage);
               phaseMaps{Cnt,iRuns} = spatialFilterGaussian(phaseMaps{Cnt,iRuns},phaseMapSmth);    
               phaseMaps{Cnt,iRuns} = phaseMaps{Cnt,iRuns}* (screenSize((iConds == 1) + 1))/2 + 61.9 ; % Translate to visual angles. Half of screenSize is used for each direction (assuming that animals eye is centered on the screen) with applying offset angle.
           else  % elevation
               a1 = mod(-angle(fTransform{iConds,iTrials}(iRuns,:,:)), 2 * pi); % angle changed from negative to positive get phase angle btw [-pi pi]
               a2 = mod(-angle(fTransform{iConds+1,iTrials}(iRuns,:,:)), 2 * pi);
               a1(isnan(a1)) = 0; a2(isnan(a2)) = 0;
               phaseMaps{Cnt,iRuns} = imrotate(squeeze((a2 - a1) / 2),rotateImage);
               phaseMaps{Cnt,iRuns} = spatialFilterGaussian(phaseMaps{Cnt,iRuns},phaseMapSmth);      
%                phaseMaps{Cnt,iRuns} = phaseMaps{Cnt,iRuns}* screenSize((iConds == 1) + 1)/2 -4.9; % offset
               phaseMaps{Cnt,iRuns} = rad2deg(phaseMaps{Cnt,iRuns})*screenSize((iConds == 1) + 1)/180 -4.9;
           end                 
           phaseMaps{Cnt,iRuns}(isnan(phaseMaps{Cnt,iRuns}(:))) = 0;

        end
        cPhaseMaps{Cnt,iTrials} = median(cat(3,phaseMaps{Cnt,:}),3);
        Cnt = Cnt+1;
    end
    cMagMaps{iTrials} = median(cat(3,magMaps{:}),3);
    cMagMaps{iTrials} =(cMagMaps{iTrials}-min(cMagMaps{iTrials}(:)))./(max(cMagMaps{iTrials}(:))- min(cMagMaps{iTrials}(:))); %normalize between 0 and 1

    % compute visual field sign maps. First compute gradients and atan - same as in the 2014 Callaway paper.
    for iRuns = 1:size(fTransform{iConds,iTrials},1)
        [dhdx, dhdy] = gradient(phaseMaps{1,iRuns}); % get gradient of horizontal phasemap
        [dvdx, dvdy] = gradient(phaseMaps{2,iRuns});
        
        graddir_hor = atan2(dhdy,dhdx); % get gradient direction angle 
        graddir_vert = atan2(dvdy,dvdx);
        vdiff = exp(-1i*graddir_hor) .* exp(1i*graddir_vert); % changed to vertical -horizontal to get blue for V1 
        
        VFS{iRuns} = sin(angle(vdiff)); %Visual field sign map
        VFS{iRuns} = spatialFilterGaussian(VFS{iRuns},phaseMapSmth);
    end
    
    cVFS{1,iTrials} = median(cat(3,VFS{:}),3);
    clear magMaps phaseMaps VFS
    
    s = figure;
    subplot(2,2,1);
    imagesc(cPhaseMaps{1,iTrials});axis image; colormap hsv; colorbar; %freezeColors;
    title(['Elevation - nTrials = ' num2str(nTrials(iTrials))]); hold on 
    caxis([-40 50])
    
    subplot(2,2,2);
    imagesc(cPhaseMaps{2,iTrials});axis image; colormap jet; colorbar; %freezeColors;
    title(['Azimuth - nTrials = ' num2str(nTrials(iTrials))]);  caxis([-60 140]); 
    
    subplot(2,2,3);
    imagesc(cMagMaps{iTrials});axis image; colorbar; colormap jet;
    title('Mean Magnitude'); 

    subplot(2,2,4);
    imagesc(spatialFilterGaussian(cVFS{1,iTrials},phaseMapSmth)); axis image;colorbar
    caxis([-0.5 0.5]) % apply colorbar threshold
    title(['VisualFieldSign - binSize = ' num2str(1) '; smth = ' num2str(phaseMapSmth)]);

    savefig(s,[img.savePath filesep Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.fig']);
    s.PaperUnits = 'inches';
    set(s, 'PaperPosition', [0 0 10 10]); %createtextbox(s); createarrow(s);
    saveas(s,[img.savePath filesep Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.png'])
    clear s
end

%% get phase and amplitude + vessel map for plotting
plotPhaseMap = spatialFilterGaussian(cVFS{1,iTrials},phaseMapSmth);
plotPhaseMap = imresize(plotPhaseMap,binSize);

plotAmpMap = spatialFilterGaussian(imresize(cMagMaps{iTrials},binSize),20); %smoothed magnitude map
plotAmpMap =(plotAmpMap-min(plotAmpMap(:)))./(max(plotAmpMap(:))- min(plotAmpMap(:))); %normalize between 0 and 1
plotAmpMap = spatialFilterGaussian(plotAmpMap,20); %smoothed magnitude map

%get vessel image.
load([img.fPath filesep 'Snapshot_1.mat'],'snap');
snap = double(imrotate(snap,rotateImage));
snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1

if size(snap,1) > size(plotPhaseMap,1) %resize if vessel map is larger
    plotPhaseMap = imresize(plotPhaseMap,size(snap,1)/size(plotPhaseMap,1));
    plotAmpMap = imresize(plotAmpMap,size(snap,1)/size(plotAmpMap,1));
end

%% plot vessel map with overlayed sign map
VFS_Outline(plotPhaseMap,cPhaseMaps{2,iTrials},cPhaseMaps{1,iTrials}, Animal, img.savePath); % plot outline with phasemaps overlaid

s = figure;
imagesc(plotPhaseMap); colormap jet; 
caxis([-0.5 0.5])
title([Animal ' - PhaseMap']); axis image
savefig(s,[img.savePath filesep Animal '_RawPhaseMap.fig']);
saveas(s,[img.savePath filesep Animal '_RawPhaseMap.jpg']);
save([img.savePath filesep 'plotPhaseMap.mat'],'plotPhaseMap');
save([img.savePath filesep 'plotAmpMap.mat'],'plotAmpMap');
save([img.savePath filesep 'cMagMaps.mat'],'cMagMaps');
save([img.savePath filesep 'cPhaseMaps.mat'],'cPhaseMaps');
save([img.savePath filesep 'cVFS.mat'],'cVFS');

s = figure;
imagesc(snap);axis image; colormap gray; freezeColors;
hold on
vfsIm = imagesc(plotPhaseMap); colormap jet; 
caxis([-0.5 0.5]);
set(vfsIm,'AlphaData',plotAmpMap);axis image; 
if img. refAlign
    f= getframe;
end 
title([Animal ' - PhaseMap + Vesselmap']) 
savefig(s,[img.savePath filesep Animal '_phaseMap.fig']);
saveas(s,[img.savePath filesep Animal '_phaseMap.png']);

%%  image alignment with Allen reference outline map (optional)
% if img.refAlign   
%     load('allenDorsalMapSM.mat')
%     snapDorsal = imshow([img.fPath filesep 'dorsal.jpg']);
% %     f = getframe;
%     effPixSize = 6./[size(f.cdata,2) size(f.cdata,1)]*1e3; % in um
% %     scaleFactor = (effPixSize./dorsalMaps.allenPixelSize).*size(dorsalMaps.edgeMap);
% % 
% %     if size(dorsalMaps.edgeMap,1) < size(f.cdata,1) %resize if vessel map is larger
% %         snapDosalScaled = imresize(f.cdata, round(scaleFactor));
% %     end
% % 
% %     d = figure; imshow(snap)
% %     f = getframe(d);
%     figure; imshow(f.cdata)   
%     d = figure; imshow(dorsalMaps.edgeMapScaled); colormap gray;
%     m = getframe;
%     savefig(d, [img.savePath filesep 'AllenDorsalEdgeMap.fig']);
%     %im1.fig is the saved original image
% 
%     % %% whole dorsal map registration to common atlas 
%     % h = figure;
%     % imagesc(snap); colormap gray;  f = getframe(h); axis image; 
%     % d = figure; imagesc(dorsalMaps.edgeMapScaled); colormap gray
%     % m = getframe(d); 
%     [f1, f_map] = frame2im(f); 
%     [m1, m_map] = frame2im(m);
%     [m1, m2, both] = imalign(f1,m1);
%     
%     imshow([img.fPath filesep 'dorsal.jpg']); a = getframe;
%     imshow(both); b=getframe;
%         [a1, f_map] = frame2im(a); 
%     [b1, m_map] = frame2im(b);
% [a, b, both2] = imalign(b1,a1);
%     figure;
%     imshow(both); axis image; colormap gray; %freezeColors;
%     title([Animal ' - PhaseMap + Vesselmap + dorsalMap'])
%     savefig(f,[img.savePath filesep Animal '_PhaseMap + Vesselmap + dorsalMap.fig']);
%     saveas(f,[img.savePath filesep Animal '_PhaseMap + Vesselmap + dorsalMap.jpg'])
%     save([img.savePath filesep Animal 'Workspace'], 'img', 'effPixSize', 'scaleFactor', 'StimDataset');
%     
% end 

%nested functions
function img = spatialFilterGaussian(img, sigma)
if sigma > 0 && (numel(img) ~=  sum(sum(isnan(img))))
    hh = fspecial('gaussian',size(img),sigma);
    hh = hh/sum(hh(:));
    img = ifft2(fft2(img).*abs(fft2(hh)));
end
end