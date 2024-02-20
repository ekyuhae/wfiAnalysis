function stimOn = findStimOnFrame(time, offset, frameTimes)
% find stimulus onset frame index 
% after stimulus timings are sorted by particular stimulus parameter, you
% input those as "time" and give some time offset to find the frameTime
% index that lies within that stimulus onset +/- offset range. 
% for 15Hz effective fps, 50ms of offset is appropriate.
% EK 2023 

frameCounter=1;
if ~isempty(time)
    for t = 1:length(time)
        timeOffset = offset ; % find frametimes within stimulus time -50ms ~ +50ms range
        stimFrames = find(frameTimes >= time(t) & frameTimes <= (time(t) + 2*timeOffset), 1, 'first');
        if ~isempty(stimFrames)
            stimOn(frameCounter) = stimFrames - 3 ;
        else %if recorded time ended earlier than the stimulus offset
            stimOn = stimOn;
        end
        frameCounter = frameCounter + 1;
    end
else
    stimOn = [];
end
end 