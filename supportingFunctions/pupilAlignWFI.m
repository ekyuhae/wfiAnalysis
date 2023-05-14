function [pupil_aligned, me_aligned] = pupilAlignWFI(photodiode_tt, stimOn, Animal, Timeline, recdate, n_ses, dlc)
% pulling preprocessed pupil/me trace and align w MC start time
%% format for pulling pupil
animal = Animal;
Date = recdate;
n_ses = num2str(n_ses);

%% finding camera start time to align pupil frames to photodiode
for h = 1:size(Timeline.hw.inputs,2)
   tp = strcmp([Timeline.hw.inputs(h).name], ['photoDiode']);
   tc = strcmp([Timeline.hw.inputs(h).name], ['eyeCameraStrobe']);
   if tp == 1
       photoIdx = h;
   end
   
   if tc == 1
       camIdx = h;
   end
end

tlPhotodiode = Timeline.rawDAQData(:,photoIdx);
eyeCam = Timeline.rawDAQData(:,camIdx);

numSamples = Timeline.rawDAQSampleCount;
fs2 = Timeline.hw.daqSampleRate;
tt_tl = 0:1/fs2:numSamples./fs2;
tt_tl = tt_tl(1,1:end-1);

threshTL = mean(tlPhotodiode);
pdt_tl = abs([diff(tlPhotodiode > threshTL)]) > 0;
bars2 = find(pdt_tl);
camStart = tt_tl(find(pdt_tl,1,'first')); % based on timeline computer 

% threshEC = mean(eyeCam);
% EC_TL = abs(diff(eyeCam > threshEC)) >0;
% frameFlips= find(EC_TL);
% frameFlips_halved = frameFlips(1:2:end);
% variable to get everything aligned to the same time-frame
brTime = stimOn;

%% figure output to verify alignment was done correctly
% figure;
% subplot(2,1,1)
% plot(tt_tl(1:end), [eyeCam(1:end)+15]); hold on;
% plot(photodiode_tt(1:end), filtered.photodiode(1:100000)./1000, 'Color', [128 128 128]./255); hold on;
% xline(photoStart, '--k'); hold on;
% title(['EyeStrobe Before Alignment']);
% xlim([0 50])
% ylim([15 23]);
% 
% subplot(2,1,2)
% plot(tt_tl(1:end)+brTime, [eyeCam(1:1000000)+15]); hold on;
% plot(photodiode_tt(1:100000), filtered.photodiode(1:100000)./1000, 'Color', [128 128 128]./255); hold on;
% title(['EyeStrobe After Alignment']);
% xline(photoStart, '--k'); hold on;
% xlim([0 50])
% ylim([15 23]);


%% pulling in pupil results from processed video
sess = n_ses;
cd(fullfile('Y:\haider\Data\EYE\ExpData\', animal, Date, sess))
if dlc == 0
    S = fullfile('Y:\haider\Data\EYE\ExpData\', animal, Date, sess, [Date '_' sess '_' animal '_eye_processed.mat']);
    load(S);
    ra = results.area;
elseif dlc == 1
    S = fullfile('Y:\haider\Data\EYE\ExpData\', animal, Date, sess, [Date '_' sess '_' animal '_eye_processed_dlc.mat']);
    load(S);
    ra = results_dlc.area;
end

% define parameters of interest
fs = 30;                %Hz (sampling over 30 ms bins)
pupil_area = ra;

%timing and alignment
pupilTime = 0:1/fs:length(pupil_area)/fs;
pupilTime = pupilTime(1, 1:end-1);
pupilTime = pupilTime + brTime; % sec

% put frame info into cohesive structure
pupil_aligned.area = pupil_area;             %30 Hz
pupil_aligned.time = pupilTime';
pupil_aligned.realidx = find(~isnan(pupil_aligned.area));
pupil_aligned.propGoodFrames = size(pupil_aligned.realidx,1)./size(pupil_aligned.area,1);

me_aligned.vector(:,1) =results.ME(~isnan(results.ME)); %skip first index which is nan
me_aligned.time = pupil_aligned.time;
save('pupilAlign.mat', 'pupil_aligned','me_aligned');

end
