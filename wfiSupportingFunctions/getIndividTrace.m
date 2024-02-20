function [dFFtrace, roi] = getIndividTrace (paddedResponseFrames, tempwindow, allData, roi, cord, sPath, flag)

% 8s window, frame rate = 30
% window = 15;
temp = reshape(paddedResponseFrames', tempwindow, [])'; %matrix to get n (number of repeats) traces 

idx = find(temp > length(allData));
if ~isempty(find(temp > length(allData)))
    if size(temp,1) > 1 
        temp(end,:) = []; % EK debugged 10/16/23
    else 
        temp(idx) = [];
    end 
end 

tmpSize = size(temp,1);
if isempty (roi) % get mean trace
    [dFFtrace(1,:),roi] = getAvgPixValROI (allData(:,:,temp(1,:)),tempwindow, [], cord, sPath,flag);
    for i = 2:tmpSize
        [dFFtrace(i,:),roi] = getAvgPixValROI (allData(:,:,temp(i,:)),tempwindow, roi, cord, sPath,flag);
    end
else
    for i = 1:size(temp,1)
        [dFFtrace(i,:),roi] = getAvgPixValROI (allData(:,:,temp(i,:)),tempwindow, roi, cord, sPath,flag);
    end
end