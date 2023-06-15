function [positionsAndTimes, paramSorted] = behSortParameter(beh_all, block, parameters, bFrameTimes, vFrameTimes, stimOn, type)
% 2 channel version. sort grating timing by corresponding azimuth location and contarasts.
% integrated version w behavioral data
% created by EK Feb23 

switch type 
    
    case 'whole'
        positionsAndTimes = zeros(length(block.trial),3);
        
        for i = 1:length(block.trial)
          
            if ~isempty(block.trial(1,i).gratingStartedTime) % if the last trial is aborted
                positionsAndTimes(i,1) = block.trial(1,i).condition.targetAzimuth;
                positionsAndTimes(i,2) = block.trial(1,i).gratingStartedTime;
                positionsAndTimes(i,3) = block.trial(1,i).condition.trialContrast;
            else
                positionsAndTimes(i,:) = [];
            end
        end
        
    case 'hit'        
        idx = find(beh_all.hits_idx==1);
        positionsAndTimes = zeros(length(idx),3);
       
    case 'miss'
        idx = find(beh_all.no_licks_idx==1);
        positionsAndTimes = zeros(length(idx),3);
        
    case 'FA'
        idx = find(beh_all.falseAlarms_idx==1);
        positionsAndTimes = zeros(length(idx),3);
        
    case 'CR'
        idx = find(beh_all.falseAlarms_rej_idx==1);
        positionsAndTimes = zeros(length(idx),3);
        
    case 'Late'
        idx = find(beh_all.late_misses_idx==1);
        positionsAndTimes = zeros(length(idx),3);       
end

if ~strcmp(type, 'whole')
     for i = 1:length(idx)
           
            if ~isempty(block.trial(1,idx(i)).gratingStartedTime)
                positionsAndTimes(i,1) = block.trial(1,idx(i)).condition.targetAzimuth;
                positionsAndTimes(i,2) = block.trial(1,idx(i)).gratingStartedTime;
                positionsAndTimes(i,3) = block.trial(1,idx(i)).condition.trialContrast;
            else
                positionsAndTimes(i,:) = [];
            end
     end
end 

if bFrameTimes (1) < vFrameTimes (1) % aligning stimulus start time to frame start time (ms)
    positionsAndTimes(:,2) = (positionsAndTimes(:,2)+stimOn)*1e3 + bFrameTimes(1);
    
else 
    positionsAndTimes(:,2) =  (positionsAndTimes(:,2)+stimOn)*1e3 + vFrameTimes(1);
end

% sort binoc - monoc with contrast levels 
bcont = parameters.targetContrast(:, find(parameters.targetAzimuth == 0,1, 'first'))'; 
bcont = sort(bcont);
mcont = parameters.targetContrast(:, find( parameters.targetAzimuth == 70,1, 'first'))'; % get contrast range 
mcont = sort(mcont);
loc = unique(parameters.targetAzimuth); % get location 
param(1,:) = repelem(loc,1,length(bcont));
param(2,:) = [bcont mcont];

for i = 1: length(param)
    idx = find(param(1,i) == positionsAndTimes(:,1) & param(2,i) == positionsAndTimes(:,3) );
    paramSorted{1,i}(1:2,1) = param(:,i); % 1st row: stim loc , 2nd row: contrast
    paramSorted{1,i}(3:2+length(idx),1) = positionsAndTimes(idx,2);
end

