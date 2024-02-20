function [photodiode_tt, photoStart, fs] = getPhotoStart_timeline(Timeline)
% code from KP
for x = 1:size(Timeline.hw.inputs,2)
    if strcmp(Timeline.hw.inputs(x).name, {'photoDiode2'}) == 1
        col = x;
    end
end

photodiode = Timeline.rawDAQData(:,col); %raw pd voltage signals
fs = Timeline.hw.daqSampleRate;
% photodiode_tt = Timeline.rawDAQTimestamps;
photodiode_tt = (0:1/fs:length(photodiode)/fs); %timestamp from 0 to last timestamp
photodiode_tt = photodiode_tt(1:end-1);
% isequal(photodiode_tt, phododiode_tt1)

thresh = nanmean(photodiode);
pdt = abs([diff(photodiode>thresh)'])>0;
bars = find(pdt);
photoStart = photodiode_tt(find(pdt,1,'first'));

end

