function [photodiode_tt, photoStart] = getPhotoStart_WFI(img, Analog)

photodiode = Analog(img.stimLine,:); %raw pd voltage signals

% figure;
% plot(1:length(Analog), Analog(2,:))

fs = 1000;

photodiode_tt = (0:1/fs:length(photodiode)/fs); %timestamp from 0 to last timestamp
photodiode_tt = photodiode_tt(1:end-1);
% isequal(photodiode_tt, phododiode_tt1)

thresh = nanmean(photodiode);
pdt = abs([diff(photodiode>thresh)'])>0;
bars = find(pdt);
photoStart = photodiode_tt(find(pdt,1));
% photoStart = photoStart(1); %first stim on
end