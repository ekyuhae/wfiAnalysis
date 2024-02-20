function [positionsAndTimes, paramSorted] = barmapSortParameter(block, parameters, bFrameTimes, vFrameTimes, stimOn)
% 2 channel version. sort bar timing corresponding to bar contrast, position
% created by EK Feb23 

positionsAndTimes = zeros(length(block.trial),3);

for i = 1:length(block.trial)
    positionsAndTimes(i,1) = block.trial(1,i).condition.position;
    positionsAndTimes(i,2) = block.trial(1,i).stimulusStartedTime;
    positionsAndTimes(i,3) = mode(block.trial(1,i).condition.colour); % reduces [r,g,b] to r but for black and white that's enough
end

if bFrameTimes (1) < vFrameTimes (1) % aligning stimulus start time to frame start time (ms)
    positionsAndTimes(:,2) = (positionsAndTimes(:,2)+stimOn)*1e3 + bFrameTimes(1);
else 
    positionsAndTimes(:,2) =  (positionsAndTimes(:,2)+stimOn)*1e3 + vFrameTimes(1);
end

param(1,:) = parameters.colour(1,:); param(2,:) = parameters.position;
paramSorted = [];
color = unique(parameters.colour);
for c = 1:length(color)
    idx = find(param(1,:)== color(c));
    paramSorted = [paramSorted param(:,idx)]; %sorted by color 
    idx = find(paramSorted(1,:)== color(c));
    paramSorted(2,idx) = sort(paramSorted(2,idx)); % sort position
end

% find corresponding stimulus onset time
rep = unique(parameters.numRepeats);
paramSorted =vertcat(paramSorted, zeros(rep,length(paramSorted)));
for i = 1:length(paramSorted)
    idx = find(paramSorted(1,i)== positionsAndTimes(:,3) & paramSorted(2,i) == positionsAndTimes(:,1));
    paramSorted(3:end, i) = positionsAndTimes(idx,2); 
end 