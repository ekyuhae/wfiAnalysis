function [avgStimResponse, avgTempTrace, binocOn, monocOn] = behSpatioTempResp_event(paramSorted, allData, bFrameTimes, sPath, fwindow,img, type, hva, reactionTriggered)
% get temporal trace and spatial activity map from the imaging during spatial detection
% task experiment.
% modified version of passive grating analysis (barmapFullcont_gratings) 
% EK Feb23
% NA 6/1/2023 - added hva as last input
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

avgStimResponse.binoc = cell(1,length(bcont));
avgStimResponse.monoc = cell(1,length(mcont));
indivTempTrace.binoc = cell(1,length(bcont));
indivTempTrace.monoc = cell(1,length(mcont));
avgTempTrace.binoc = zeros(length(bcont), fr*3+1);
avgTempTrace.monoc = zeros(length(mcont), fr*3+1);

for cont = 1:length(bcont) 
    %% spatial analysis
        binocOn{cont} = findStimOnFrame(binoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(binocOn{cont})
            for i = 1:length(binocOn{cont})
                binocFrames = cat(2,binocFrames, binocOn{cont}(i):binocOn{cont}(i)+fwindow);
                binocTemp= cat(2,binocTemp, binocOn{cont}(i)-fr:binocOn{cont}(i)+fr*2); % 6s window - 2s prestim, 4s poststim
            end
            binocIndFrame{cont} = allData(:,:,binocFrames);
            avgStimResponse.binoc{cont} = nanmean(binocIndFrame{cont}, 3);
        end
        btraceLength = length(binocTemp) / length(binocOn{cont});

        monocOn{cont} = findStimOnFrame(monoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(monocOn{cont})
            for i = 1:length(monocOn{cont})
                monocFrames = cat(2,monocFrames, monocOn{cont}(i):monocOn{cont}(i)+fwindow);
                monocTemp= cat(2,monocTemp, monocOn{cont}(i)-fr:monocOn{cont}(i)+fr*2);
            end
            monocIndFrame{cont} = allData(:,:,monocFrames);
            avgStimResponse.monoc{cont} = nanmean(monocIndFrame{cont}, 3);
        end
        mtraceLength= length(monocTemp) / length(monocOn{cont});
        %% temporal analysis - requires alignment coordinates and ROI pixels 
        if cont ==1
            if isempty(hva)
                if ~reactionTriggered
                    cfile =[sPath '\..' filesep num2str(1) filesep 'analyzed_active_1.mat'];
                else 
                    cfile =[sPath '\..\..' filesep num2str(1) filesep 'analyzed_active_1.mat'];
                end 
            else
                cfile = [sPath '\..\..'  filesep num2str(1) filesep hva filesep 'analyzed_active_1_' hva '.mat'];
            end
            if exist(cfile)
                load(cfile,'roi', 'cord')
                disp('roi data has been loaded')
            end
        end
        if ~isempty(hva)
            sPath= fullfile(sPath, '..');
        end 
        if ~isempty(binocTemp)
%             [avgTempTrace.binoc(cont,:),~] = getAvgTemporalResponse(binocTemp, btraceLength, allData, roi.binoc,cord, sPath); %AR Changed
            indivTempTrace.binoc{cont} = getIndividTrace(binocTemp,btraceLength, allData, roi.binoc, cord, sPath, flag);
            avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
        end
        if ~isempty(monocTemp)
%             [avgTempTrace.monoc(cont,:),~] = getAvgTemporalResponse(monocTemp, mtraceLength, allData, roi.monoc,cord, sPath); %AR Changed
            indivTempTrace.monoc{cont} = getIndividTrace(monocTemp,mtraceLength, allData, roi.monoc, cord, sPath, flag);
            avgTempTrace.monoc(cont,:) = mean(indivTempTrace.monoc{cont},1);
        end
        binocFrames = [];
        binocTemp = [];
        monocFrames = [];
        monocTemp = [];
        if ~isempty(hva)
            sPath = sPath(1:end-2);
        end 
end 

cd(sPath)

%if isempty (hva)
plotAvgMap_beh(contrasts, avgStimResponse,cord, sPath, type, hva, reactionTriggered) % contrast range for binoc and monoc is different 5/11 fix ity
% else 
%    plotAvgMap_beh(contrasts, avgStimResponse,cord, fullfile(sPath, '..'), type, hva, reactionTriggered) % contrast range for binoc and monoc is different 5/11 fix ity 
% end 

plotTempTrace_beh(avgTempTrace, contrasts, fr, type, hva) 

if isempty(hva)
    if ~reactionTriggered
        save(['analyzed_active_' num2str(sPath(end), '%04i') '_' type], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
    else
        save(['analyzed_active_' num2str(sPath(end-6), '%04i') '_' type], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
    end
else
    if ~reactionTriggered
        sessDir = fileparts(pwd);
        save(['analyzed_active_' num2str(sessDir(end)) '_' hva '_' type], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
    else
        sessDir = fileparts(pwd);
        save(['analyzed_active_' num2str(sessDir(end-6)) '_' hva '_' type], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
    end
end

%% 
end