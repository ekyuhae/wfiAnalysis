%% WFI eyecam preprocessed data 
% get pupil size and ME structure
% date = '2023-02-01';
% sessNr ='4'
fPath = 'Y:\haider\Data\EYE\ExpData\M230131_1\2023-03-01\4';
addpath (genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Pupil_Analysis'));
addpath (genpath(fPath));
results = getResultsStructure(['2023-03-01_4_M230131_1_eyeDLC_resnet50_PupilDLCZooModelJun22shuffle1_85000_bx.csv']);

save([fPath filesep '2023-03-01_4_M230131_1_eye_processed.mat'], 'results');

figure; 
subplot(2,1,1)
plot (results.ME)
xlabel('frame'); ylabel('area(pixels)')
title(['ME trace: ' fPath(28:end)])
subplot(2,1,2)
plot (results.area)
xlabel('frame'); ylabel('area(pixels)')
title(['Pupil area trace: ' fPath(28:end)])
