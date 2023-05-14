function [bdataCorr, vdataCorr] = motionCorrect(bData, vData)
% perform WFI motion correction 
% EK Feb23

    bData = single(squeeze(bData));
    blueRef = fft2(median(bData,3)); %blue reference for motion correction
    
    vData = single(squeeze(vData));
    violetRef = fft2(median(vData,3)); %violet reference for motion correction
        
%perform motion correction for both channels
for iFrames = 1:size(bData,3)
    [~, temp] = dftregistration(blueRef, fft2(bData(:, :, iFrames)), 10);
    bdataCorr(:, :, iFrames) = abs(ifft2(temp));
    
    [~, temp] = dftregistration(violetRef, fft2(vData(:, :, iFrames)), 10);
    vdataCorr(:, :, iFrames) = abs(ifft2(temp));
end
end 