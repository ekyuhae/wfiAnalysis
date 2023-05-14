function [photodiode_tt, photoStart] = getPhotoStartFixed_WFI(img, Analog)

photodiode = Analog(img.stimLine,:); %raw pd voltage signals

% figure;
% plot(1:length(Analog), Analog(5,:))

fs = 1000;

photodiode_tt = (0:1/fs:length(photodiode)/fs); %timestamp from 0 to last timestamp
photodiode_tt = photodiode_tt(1:end-1);
% isequal(photodiode_tt, phododiode_tt1)
ma = max(photodiode); mi = min(photodiode);
if photodiode(1) > 4000 % starts at high signal
    thresh = 0.5*(ma+mi);
elseif photodiode(1) - mi < 500
    thresh = nanmean(photodiode);
elseif photodiode(1) < nanmean(photodiode) % starts at middle
    thresh = 0.5*(photodiode(1) + mi);
end 

pdt = abs([diff(photodiode>thresh)'])>0;
bars = find(pdt);
photoStart = photodiode_tt(find(pdt,1));
% photoStart = photoStart(2); %first stim on
end