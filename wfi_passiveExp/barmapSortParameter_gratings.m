function [positionsAndTimes, paramSorted] = barmapSortParameter_gratings(block, parameters, bFrameTimes, vFrameTimes, stimOn)
% 2 channel version. sort grating timing by corresponding azimuth location.
% full contrast ver
% created by EK Feb23, NA July23 

if length(stimOn) == 2
    stimOn = stimOn(2);
end

positionsAndTimes = zeros(length(block.trial),3);

for i = 1:length(block.trial)
    positionsAndTimes(i,1) = block.trial(1,i).condition.targetAzimuth;
    positionsAndTimes(i,2) = block.trial(1,i).gratingStartedTime;
    positionsAndTimes(i,3) = block.trial(1,i).condition.trialContrast; % reduces [r,g,b] to r but for black and white that's enough
end

if bFrameTimes (1) < vFrameTimes (1) % aligning stimulus start time to frame start time (ms)
    positionsAndTimes(:,2) = (positionsAndTimes(:,2)+stimOn)*1e3 + bFrameTimes(1);
    
else 
    positionsAndTimes(:,2) =  (positionsAndTimes(:,2)+stimOn)*1e3 + vFrameTimes(1);
end

cont = unique(parameters.targetContrast(:,1))'; % get contrast range 
loc = unique(parameters.targetAzimuth); % get location 
param = [repelem(loc,1,length(cont)); repmat(cont,1,length(loc))];

for i = 1: length(param)
    idx = find(param(1,i) == positionsAndTimes(:,1) & param(2,i) == positionsAndTimes(:,3) );
    paramSorted{1,i}(1:2,1) = param(:,i);
    paramSorted{1,i}(3:2+length(idx),1) = positionsAndTimes(idx,2);
end

