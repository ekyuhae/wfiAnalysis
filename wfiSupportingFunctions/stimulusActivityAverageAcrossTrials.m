
function AzimuthAverages = trialAverageResponse(trialFrame)
%stimulusActivityAverageAcrossTrials -- take trialFrame from
%PhaseMapAnalysis pipeline and average each trial (3-D array) along each sweep
%direction in the trialFrame cell array.
% output | compositeAzimuthAverages (2-D cell array): trial averaged
% frames.

sweepDirectionCount = size(trialFrame, 1);
trialCount = size(trialFrame, 2);

if ~(isempty(sweepDirectionCount) && isempty(trialCount))
    for sweepDirectionNumber = 1:sweepDirectionCount
    azimuthAverage = [];
        for trial = 1:trialCount
            azimuthAverage = cat(4,azimuthAverage, trialFrame{sweepDirectionNumber, trial});
        end
        AzimuthAverages{sweepDirectionNumber, 1} = mean(azimuthAverage, 4);
    end
end
end

