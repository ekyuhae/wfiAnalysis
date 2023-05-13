function [avgStimResponse, avgTempTrace, binocOn, monocOn] = behSpatioTempResp_event(paramSorted, allData, bFrameTimes, sPath, fwindow,img, type)
% get temporal trace and spatial activity map from the imaging during spatial detection
% task experiment.
% modified version of passive grating analysis (barmapFullcont_gratings) 
% EK Feb23
%% 

for i = 1:length(paramSorted)
    contrasts(i) = paramSorted{1,i}(2);
end 
idx = find(contrasts == 0, 1, 'last');
bcont = contrasts(1:idx-1);
mcont = contrasts(idx:end);

binoc = paramSorted(1:idx-1);
monoc = paramSorted(idx:end);

%% WFI barmapping analysis

binocOn = [];
monocOn = [];
binocFrames = [];
binocTemp = [];
monocFrames = []; 
monocTemp = [];
fr = img.sRate;

for cont = 1:length(bcont) 
    %% spatial analysis
        binocOn{cont} = findStimOnFrame(binoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(binocOn{cont})
            for i = 1:length(binocOn{cont})
                binocFrames = cat(2,binocFrames, binocOn{cont}(i):binocOn{cont}(i)+fwindow);
                binocTemp= cat(2,binocTemp, binocOn{cont}(i)-fr:binocOn{cont}(i)+fr*1.5); % 5s window
            end
            binocIndFrame{cont} = allData(:,:,binocFrames);
            avgStimResponse.binoc{cont} = nanmean(binocIndFrame{cont}, 3);
        end
        btraceLength = length(binocTemp) / length(binocOn{cont});

        monocOn{cont} = findStimOnFrame(monoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(monocOn{cont})
            for i = 1:length(monocOn{cont})
                monocFrames = cat(2,monocFrames, monocOn{cont}(i):monocOn{cont}(i)+fwindow);
                monocTemp= cat(2,monocTemp, monocOn{cont}(i)-fr:monocOn{cont}(i)+fr*1.5);
            end
            monocIndFrame{cont} = allData(:,:,monocFrames);
            avgStimResponse.monoc{cont} = nanmean(monocIndFrame{cont}, 3);
        end
        mtraceLength= length(monocTemp) / length(monocOn{cont});
        %% temporal analysis
      
        cfile =[sPath(1:end-1) filesep num2str(1) filesep 'analyzed_active_1.mat'];
        if exist(cfile)
            load(cfile,'roi', 'cord')
            disp('roi data has been loaded')
        end

        [avgTempTrace.binoc(cont,:),~] = getAvgTemporalResponse(binocTemp, btraceLength, allData, roi.binoc,cord, sPath); %AR Changed
        [avgTempTrace.monoc(cont,:),~] = getAvgTemporalResponse(monocTemp, mtraceLength, allData, roi.monoc,cord, sPath); %AR Changed

        indivTempTrace.binoc{cont} = getIndividTrace(binocTemp,btraceLength, allData, roi.binoc, cord, sPath);
        indivTempTrace.monoc{cont} = getIndividTrace(monocTemp,mtraceLength, allData, roi.monoc, cord, sPath);

        binocFrames = [];
        binocTemp = [];
        monocFrames = [];
        monocTemp = [];
end 

cd(sPath)

plotAvgMap_beh(contrasts, avgStimResponse,cord, sPath) % contrast range for binoc and monoc is different 5/11 fix ity
plotTempTrace_beh(avgTempTrace, contrasts, fr) 

save(['analyzed_active_' num2str(sPath(end), '%04i') '_' type], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')


%% 
end