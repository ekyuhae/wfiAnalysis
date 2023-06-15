function [avgStimResponse, avgTempTrace, binocOn, monocOn] = behSpatioTempResp(paramSorted, allData, bFrameTimes, sPath, fwindow,img, sessNr, expID_ret, hva, flag)
% get temporal trace and spatial activity map from the imaging during spatial detection
% task experiment.
% modified version of passive grating analysis (barmapFullcont_gratings) 
% EK Feb23
% NA 6/1/2023 - added in hva input and change in plots based on this
%% 

for i = 1:length(paramSorted)
    contrasts(i) = paramSorted{1,i}(2);
end 
idx = find(contrasts == 0, 1, 'last');
bcont = contrasts(1:idx-1);
mcont = contrasts(idx:end);

binocc = paramSorted(1:idx-1);
monocc = paramSorted(idx:end);

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
        binocOn{cont} = findStimOnFrame(binocc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(binocOn{cont})
            for i = 1:length(binocOn{cont})
                binocFrames = cat(2,binocFrames, binocOn{cont}(i):binocOn{cont}(i)+fwindow);
                binocTemp= cat(2,binocTemp, binocOn{cont}(i)-fr:binocOn{cont}(i)+fr*2); % 5s window
            end
            binocIndFrame{cont} = allData(:,:,binocFrames);
            avgStimResponse.binoc{cont} = nanmean(binocIndFrame{cont}, 3);
        end
        btraceLength = length(binocTemp) / length(binocOn{cont});

        monocOn{cont} = findStimOnFrame(monocc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(monocOn{cont})
            for i = 1:length(monocOn{cont})
                monocFrames = cat(2,monocFrames, monocOn{cont}(i):monocOn{cont}(i)+fwindow);
                monocTemp= cat(2,monocTemp, monocOn{cont}(i)-fr:monocOn{cont}(i)+fr*2);
            end
            monocIndFrame{cont} = allData(:,:,monocFrames);
            avgStimResponse.monoc{cont} = nanmean(monocIndFrame{cont}, 3);
        end
        mtraceLength= length(monocTemp) / length(monocOn{cont});
        %% temporal analysis
        if cont == 1 && sessNr ==1
            if ~isempty(hva)
                sPath = fullfile(sPath, '..');
            end
            if exist(fullfile(sPath, 'coordinates.mat'), 'file') == 2
                load(fullfile(sPath, 'coordinates.mat'));
                disp(['coordinates.mat has been loaded']);
            else
                cord = registerRetMap(sPath, 2, expID_ret, flag); % for behavior exp (2), register retinotopy
            end
            if ~isempty(binocTemp)
                [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData, [], cord, sPath, 'binoc');
                avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
            end
            if ~isempty(monocTemp)
                [indivTempTrace.monoc{cont}, roi.monoc] = getIndividTrace(monocTemp,mtraceLength, allData, [], cord, sPath, 'monoc');
                avgTempTrace.monoc(cont,:) = mean(indivTempTrace.monoc{cont},1);
            end
        else
            if ~isempty(hva)
                sPath= fullfile(sPath, '..');
                cfile =[sPath '\..\' num2str(1) filesep hva filesep 'analyzed_active_1_' hva '.mat'];
            else
                cfile =[sPath(1:end-1) filesep num2str(1) filesep 'analyzed_active_1.mat'];
            end
            if exist(cfile) % for session second session and so on, automatically load roi that was chosen before and get avg roi trace
                load(cfile,'roi', 'cord')
                % disp('ROI has been loaded')
                if ~isempty(binocTemp)
                    if ~isfield(roi, 'binoc')
                        [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,[], cord, sPath, 'binoc');
                    else
                        [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,roi.binoc, cord, sPath,'binoc');
                    end
                    avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
                end
                if ~isempty(monocTemp)
                    if ~isfield(roi, 'monoc')
                        [indivTempTrace.monoc{cont}, roi.monoc] = getIndividTrace(monocTemp,mtraceLength, allData,[], cord, sPath, 'monoc');
                    else
                        [indivTempTrace.monoc{cont}, roi.monoc] = getIndividTrace(monocTemp,mtraceLength, allData, roi.monoc, cord, sPath, 'monoc');
                    end
                    avgTempTrace.monoc(cont,:) = mean(indivTempTrace.monoc{cont},1);
                end
            else
                if ~isempty(binocTemp)
                    if ~isfield(roi, 'binoc')
                        [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,[], cord, sPath, 'binoc');
                    else
                        [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,roi.binoc, cord, sPath, 'binoc');
                    end
                    avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
                    if ~isempty(monocTemp)
                        if ~isfield(roi, 'monoc')
                            [indivTempTrace.binoc{cont}, roi.monoc] = getIndividTrace(monocTemp,mtraceLength, allData,[], cord, sPath, 'monoc');
                        else
                            [indivTempTrace.monoc{cont}, roi.monoc] = getIndividTrace(monocTemp,mtraceLength, allData, roi.monoc, cord, sPath, 'monoc');
                        end
                        avgTempTrace.monoc(cont,:) = mean(indivTempTrace.monoc{cont},1);
                    end
                end           
            end
            %  if ~isempty(hva)
            %     sPath= fullfile(sPath, '..');
            % end
          
        %     if ~isempty(binocTemp)
        %         [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData, roi.binoc, cord, sPath,hva);
        %         avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
        %     end
        % end
        % if ~isempty(roi.monoc)
        %     if ~isempty(monocTemp)
        %         [indivTempTrace.monoc{cont}, roi.monoc] = getIndividTrace(monocTemp,mtraceLength, allData, [], cord, sPath,hva);
        %         avgTempTrace.monoc(cont,:) = mean(indivTempTrace.monoc{cont},1);
        %     end
        % else
        %     if ~isempty(monocTemp)
        %         [indivTempTrace.monoc{cont}, roi.monoc] = getIndividTrace(monocTemp,mtraceLength, allData, roi.monoc, cord, sPath, hva);
        %         avgTempTrace.monoc(cont,:) = mean(indivTempTrace.monoc{cont},1);
        %     end
        % end
        binocFrames = [];
        binocTemp = [];
        monocFrames = [];
        monocTemp = [];
        if ~isempty(hva)
            sPath = sPath(1:end-2);
        end 
end 
% if ~isempty(hva)
% sPath = [sPath(1:end-2)];
% end 

cd(sPath)
plotAvgMap_beh(contrasts, avgStimResponse,cord, sPath, type, hva, reactionTriggered) % contrast range for binoc and monoc is different 5/11 fix ity

% if isempty (hva)
% plotAvgMap_beh(contrasts, avgStimResponse,cord, sPath, [], hva) % contrast range for binoc and monoc is different 5/11 fix ity
% else
% 
%    plotAvgMap_beh(contrasts, avgStimResponse,cord, fullfile(sPath, '..'), [], hva) % contrast range for binoc and monoc is different 5/11 fix ity 
% end 

plotTempTrace_beh(avgTempTrace, contrasts, fr,[], hva) 

if isempty(hva)
save(['analyzed_active_' num2str(sessNr) '%04i'], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
else
    save(['analyzed_active_' num2str(sessNr) '_' hva], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
end

%% 
end