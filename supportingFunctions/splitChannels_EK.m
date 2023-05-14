function [blueData,blueTimes,hemoData,hemoTimes,frameTimes,falseAlign,sRate] = splitChannels_EK(opts,trialNr,fileType,Analog)
% Separating blue/violet frames from widefield imging data. Requires trigger channels for
% blue/violet LEDs and some dark frames at the end. Does not use
% stimulus information here and only returns avg channel frames. 

if ~exist('fileType', 'var') || isempty(fileType)
    fileType = '.dat';
end
falseAlign = false;

%% load data and check channel identity

%08/07/22 EK modified 
load([opts.fPath filesep 'frameTimes_' num2str(trialNr, '%04i') '.mat'], 'imgSize', 'frameTimes'); %get data size
% % AbsTime = cell2mat({allFrameTimes.AbsTime}.'); 
% % frameTimes = datenum(AbsTime) * 86400*1e3; %convert to miliseconds
frameTimes = frameTimes * 86400*1e3; 
cFile = [opts.fPath filesep opts.fName '_' num2str(trialNr, '%04i') fileType]; %current file to be read
[~, data] = loadRawData(cFile,'Frames',[], imgSize); %load video data

%reshape data to compute mean frame intensities
dSize = size(data);
data = reshape(data,[],dSize(end)); % combine x/y pixel values as row vector 
temp = zscore(mean(single(data))); % get zscore of each pixel values from mean intensities
data = squeeze(reshape(data,dSize));

% temp(end:-1:end-4) =[]; %remove extra five dark frames
% bFrame = find(temp < min(temp)*.75); %index for black frames --> violet frames EK
% if bFrame(1) == 1 %if first frame is dark, remove initial frames from data until LEDs are on
%     %remove initial dark frames
%     cIdx = find(diff(bFrame) > 1, 1); 
%     data(:,:,1:cIdx) = [];
%     dSize = size(data);
%     temp(1:cIdx) = [];
%     frameTimes(1:cIdx) = [];
% end

%determine imaging rate - either given as input or determined from data
if isfield(opts,'sRate')
    sRate = opts.sRate; %EK changed variable name
else
    sRate = 1000/(mean(diff(frameTimes))*2);
end

% check if pre- and poststim are given. use all frames if not.
% if ~isfield(opts,'preStim') || ~isfield(opts,'postStim')
%     opts.preStim = 0;
%     opts.postStim = inf;
% else
%     opts.preStim = ceil(opts.preStim * sRate); % EK number of frames 
%     opts.postStim = ceil(opts.postStim * sRate);
% end
% fs = 1000; % daq sampling rate
% actualAnalTimestamps = ceil(fs*(imgSize(4)-5)/sRate); % timepoint from analog trace that corresponds to the last frame exposure signal (should be close to falling edge)
% a = 1:length(Analog); a1 = frameTimes(1)+a; 
% lag = a1(end) - frameTimes(end-5);
% disp(['lag between the last timestamp of aligned analog trace and the actual frametime of last frame: ' num2str(lag) 'ms']) %excluding extra frames
% lag2 = length(Analog) - lag; 
% disp(['last frame timepoints from analog signal - last frametime = ' num2str(actualAnalTimestamps - lag2)])
% 
% if temp(1) > 1 %started w blue frame
%     thresh_bIdx = find(abs(diff(Analog(2,:))) > mean(Analog(2,:))); %falling edge
%     thresh_vIdx = find(diff(Analog(3,:)) > mean(Analog(3,:)));
% else 
%     thresh_bIdx = find(diff(Analog(2,:)) > mean(Analog(2,:))); 
%     thresh_vIdx = find(-diff(Analog(3,:)) > mean(Analog(3,:)));
% end 
% temp(end:-1:end-4) = [];
% find(thresh_bIdx, length(temp))
% thresh_bIdx(round(diff(thresh_bIdx))<= (sRate/2)) = [];
% thresh_vIdx(round(diff(thresh_vIdx))<= (sRate/2)) = [];
% cutBlue = thresh_bIdx(1:length(find(temp > 0)));
% cutviol = thresh_vIdx(1:length(find(temp < 0)));

if any(~isnan(opts.trigLine)) || any(opts.trigLine > size(Analog,1))
    
    trace1 = Analog(opts.trigLine,:); %blue and violet light trigger channels
%     trace2 = zscore(trace1(1,end:-1:1) - trace1(2,end:-1:1)); %invert and subtract to check color of last frame % subtract triglines and get z scores from that
   % EK 12/23/22
    trace = zscore(trace1(1,:) - trace1(2,:)); %blue-violet
    trace(round(diff(trace)) ~= 0) = 0; %don't use triggers that are only 1ms long
  
%     bIdx = find(diff(trace > 1) < 0, length(temp)/2); % blue's falling edge
%     vIdx = find(diff(trace < -1) < 0, length(temp)/2); % violet's falling edge
%     if bIdx(end) > vIdx(end)
%         trace1 = trace1 (:, 1:bIdx(end));
%     else
%         trace1 = trace1 (:, 1:vIdx(end));
%     end
% 
%     trace1 = zscore(trace1(1,end:-1:1) - trace1(2,end:-1:1));
%     lastBlue = find(trace1 > 1, 1); % find 1st trigger change point starting from the last timestamps, i.e. last blue frame point
%     lastHemo = find(trace1 <-1, 1); 
%     
%     blueLast = lastBlue < lastHemo;
%     if isempty(lastBlue) || isempty(lastHemo)
%         warning(['Failed to find trigger signals. lastBlue: ' num2str(lastBlue) '; lastHemo: ' num2str(lastHemo) '; trialNr: ' num2str(trialNr)])
%     end
   
    bFrame = find(temp < min(temp)*.75); %index for first black frame (assuming the last frame is really a dark frame) EK
%     bFrame(bFrame < round(size(temp,2) / 2)) = []; %make sure black frame is in the second half of recording.
   
%     if isempty(bFrame); bFrame = length(temp); end %if there are no black frames, just use the last one
%     bFrame = bFrame(1);

    blueInd = true(1,length(temp));
    blueInd(bFrame) = false;
       
%     blueInd = false(1,length(temp));
%     if blueLast %last frame before black is blue
%         if rem(bFrame,2) == 0 %blue frames (bFrame - 1) have uneven numbers
%             blueInd(1:2:end) = true;
%         else            
%             blueInd(2:2:end) = true;
%         end
% %         lastFrame = size(trace1,2) - lastBlue; %index for end of last frame (sample timepoints)
%         lastFrame = size(trace1,2);
%     else %last frame before black is violet
%         if rem(bFrame,2) == 0 %blue frames (bFrame - 2) have even numbers
%             blueInd(1:2:end) = true;
%         else
%             blueInd(2:2:end) = true;
%         end
%         lastFrame = size(trace1,2); %index for end of last frame
%     end
    
%     %get number of rejected frames before data was saved
%     nrFrames = bFrame - 1; %number of exposed frames in the data
%     nrTriggers = size(find(diff(Analog(opts.trigLine(1),:))> 2500),2) + size(find(diff(Analog(opts.trigLine(2),:))> 2500),2); %nr of triggers
%     removedFrames = nrTriggers - nrFrames;
%     save([opts.fPath filesep 'frameTimes_' num2str(trialNr, '%04i') '.mat'], 'imgSize', 'frameTimes', 'removedFrames'); %get data size

    %realign frameTimes based on time of last non-dark frame
%     frameTimes = (frameTimes - frameTimes(bFrame - 1)) + lastFrame;
%     blueInd = blueInd(frameTimes < size(trace1,2));
%     blueInd(bFrame - 1:end) = []; %exclude black and last non-black frame
    
    blueTimes = frameTimes(blueInd);
    hemoTimes = frameTimes(~blueInd);
    
    blueData = data(:,:,blueInd);
    hemoData = data(:,:,~blueInd);

   chanDiff = size(blueData,3) - size(hemoData,3);
    if chanDiff < 0 %less blue frames, cut hemo frame
        hemoData = hemoData(:,:,1:end+chanDiff);
        hemoTimes = hemoTimes(1:end+chanDiff);
    elseif chanDiff > 0 %less hemo frames, cut blue frame
        blueData = blueData(:,:,1:end-chanDiff);
        blueTimes = blueTimes(1:end-chanDiff);
    end
    
    if opts.plotChans
    figure %show result
    subplot(1,2,1); colormap gray
    imagesc(mean(blueData,3)); axis image
    title(['Blue frame average - Trial ' num2str(trialNr)])
    subplot(1,2,2);
    imagesc(mean(hemoData,3)); axis image
    title(['Hemo frame average - Trial ' num2str(trialNr)])
    drawnow;
    end
    
end