
function azimuthAverages = trialAveragedResponse(trialFrame)
%stimulusActivityAverageAcrossTrials -- take trialFrame from
%PhaseMapAnalysis pipeline and average each trial (3-D array) along each sweep
%direction in the trialFrame cell array.
% output | compositeAzimuthAverages (2-D cell array): trial averaged
% frames.

sweepDirectionCount = size(trialFrame, 1);
trialCount = size(trialFrame, 2);
azimuthAverages = cell(sweepDirectionCount, 1);

if ~(isempty(sweepDirectionCount) && isempty(trialCount))
    % by row of trialFrame, which is azimuth sweep or elevation sweep
    % trialFrame can have up to four (back and forth for azimuth and
    % elevation)
    for sweepDirectionNumber = 1:sweepDirectionCount
        azimuthAverage = [];
        for trial = 1:trialCount
            azimuthAverage = cat(4,azimuthAverage, trialFrame{sweepDirectionNumber, trial});
        end
        azimuthAverages{sweepDirectionNumber, 1} = mean(azimuthAverage, 4);
    end
end
end

