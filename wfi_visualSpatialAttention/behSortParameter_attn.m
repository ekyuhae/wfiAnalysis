function [bHitTimes, mHitTimes] = behSortParameter_attn(trTable, block, bFrameTimes, vFrameTimes, stimOn)
% Output: bHItTimes & mHitTimes, stimulus onset during hit trials and aligned to imaging
% FrameTimes 
% created by EK Nov23

% find hit trials in each block from trTable
bTrNr = max(trTable.trNumInBlock (trTable.boolBinoc));
mTrNr = max(trTable.trNumInBlock (trTable.boolMonoc));

idxb = cell(1,bTrNr); idxm = cell(1,mTrNr);

for i = 1: bTrNr
idxb{i} =  find(trTable.hitTrNumInBlock== i & trTable.boolBinoc);
end  

for i = 1: mTrNr
idxm{i} =  find(trTable.hitTrNumInBlock== i & trTable.boolMonoc);
end  

% find grating started times from the hit trial alinged indices for each block
for i = 1:bTrNr
    if ~isempty(idxb{i})
        for ii = 1:length(idxb{i})
            bHitTimes{i}(ii) = block.trial(1,idxb{i}(ii)).gratingStartedTime;
        end
    end
end
for i = 1:mTrNr
    if ~isempty(idxm{i})
        for ii = 1:length(idxm{i})
            mHitTimes{i}(ii) = block.trial(1,idxm{i}(ii)).gratingStartedTime;
        end
    end
end

if bFrameTimes (1) < vFrameTimes (1) % aligning stimulus start time to frame start time (ms)
    bHitTimes = cellfun(@(x) (x(:) +stimOn)*1e3 + bFrameTimes(1), bHitTimes, 'UniformOutput', false);
    mHitTimes = cellfun(@(x) (x(:) +stimOn)*1e3 + bFrameTimes(1), mHitTimes, 'UniformOutput', false);
else
    bHitTimes = cellfun(@(x) (x(:) +stimOn)*1e3 + vFrameTimes(1), bHitTimes, 'UniformOutput', false);
    mHitTimes = cellfun(@(x) (x(:) +stimOn)*1e3 + vFrameTimes(1), mHitTimes, 'UniformOutput', false);
end

% % sort binoc - monoc with contrast levels
% %bcont = parameters.targetContrast(:, find(parameters.targetAzimuth == 0,1, 'first'))';
% bcont = parameters.targetContrast(:, find(parameters.targetAzimuth == 0,1, 'first'))'; %CHANGE during phase 1  and 2 (when stimulus moving)
% bcont = sort(bcont);
% mcont = parameters.targetContrast(:, find( parameters.targetAzimuth == 70,1, 'first'))'; % get contrast range
% mcont = sort(mcont);
% loc = unique(parameters.targetAzimuth); % get location
% param(1,:) = repelem(loc,1,length(bcont));
% param(2,:) = [bcont mcont]; %change if only location is 70
% 
% for i = 1: length(param)
%     idx = find(param(1,i) == positionsAndTimes(:,1) & param(2,i) == positionsAndTimes(:,3) ); %specific location and contrast
%     paramSorted{1,i}(1:2,1) = param(:,i); % 1st row: stim loc , 2nd row: contrast
%     paramSorted{1,i}(3:2+length(idx),1) = positionsAndTimes(idx,2);
% end
% end

