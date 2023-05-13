function [avgStimResponse, avgTempTrace, binocOn, monocOn] = barmapMultiCont_gratings(paramSorted, allData, bFrameTimes, sPath, window,img, sessNr, expID_ret)
% bar mapping full contrast analysis. updated version of WFIbarmappingAnalysis2
% plot activity maps corresponding to different bar locations
% version2 of barmapFullcont_gratings, can apply for multi contrast exp too
% EK Feb23
%% 

% addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Arvind_Ramesh\AR_function'))
% addpath(genpath('\\ad.gatech.edu\bme\labs\haider\Data\analyzedData\EKK\WFI\Codes\Arvind_Ramesh\AR_function'))

for i = 1:length(paramSorted)
    contrasts(i) = paramSorted{i}(2);
end 
contrasts = unique(contrasts);

binoc = paramSorted(1:length(contrasts));
monoc = paramSorted(length(contrasts)+1:end);

%% WFI barmapping analysis

binocOn = [];
monocOn = [];
binocFrames = [];
binocTemp = [];
monocFrames = []; 
monocTemp = [];
fr = img.sRate*2;

for cont = 1:length(contrasts) 
    %% spatial analysis
        binocOn{cont} = findStimOnFrame(binoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(binocOn{cont})
            for i = 1:length(binocOn{cont})
                binocFrames = cat(2,binocFrames, binocOn{cont}(i):binocOn{cont}(i)+window);
                binocTemp= cat(2,binocTemp, binocOn{cont}(i)-fr:binocOn{cont}(i)+fr);
            end
            avgStimResponse.binoc{cont} = nanmean(allData(:, :, binocFrames), 3);
        end
        
        monocOn{cont} = findStimOnFrame(monoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(monocOn{cont})
            for i = 1:length(monocOn{cont})
                monocFrames = cat(2,monocFrames, monocOn{cont}(i):monocOn{cont}(i)+window);
                monocTemp= cat(2,monocTemp, monocOn{cont}(i)-fr:monocOn{cont}(i)+fr);
            end
            avgStimResponse.monoc{cont} = nanmean(allData(:, :, monocFrames), 3);
        end
        
        %% temporal analysis
        if cont == 1 && sessNr ==1
            if exist(fullfile(sPath, 'coordinates.mat'), 'file') == 2
                load(fullfile(sPath, 'coordinates.mat'));
                disp(['coordinates.mat has been loaded']);
            else
                cord = registerRetMap(sPath, 0, expID_ret); % for grating exp
            end
            [avgTempTrace.binoc(cont,:), roi.binoc] = getAvgTemporalResponse(binocTemp, img, allData,[],cord, sPath); %AR Changed (added sPath)
            [avgTempTrace.monoc(cont,:), roi.monoc] = getAvgTemporalResponse(monocTemp, img, allData,[],cord, sPath); %AR Changed
        else
            cfile =[sPath(1:end-1) filesep num2str(1) filesep 'analyzed_passive_1.mat'];
               if exist(cfile)
                   load(cfile,'roi', 'cord')
               end 
            
            [avgTempTrace.binoc(cont,:),~] = getAvgTemporalResponse(binocTemp, img, allData, roi.binoc,cord, sPath); %AR Changed
            [avgTempTrace.monoc(cont,:),~] = getAvgTemporalResponse(monocTemp, img, allData, roi.monoc,cord, sPath); %AR Changed
        end 
        

        binocFrames = [];
        binocTemp = [];
        monocFrames = [];
        monocTemp = [];
end 

cd(sPath)
% for i = 1: length(contrasts)
% avgStimResponse.binoc{1,i}{1,1}.cdata= imresize(avgStimResponse.binoc{1,i}{1,1}.cdata, [size(allData,1) size(allData,2)]);
% avgStimResponse.monoc{1,i}{1,1}.cdata= imresize(avgStimResponse.monoc{1,i}{1,1}.cdata, [size(allData,1) size(allData,2)]);
% end
plotAvgMap_gratings(contrasts, avgStimResponse,cord, sPath)
plotTempTrace_grating(avgTempTrace, contrasts, fr) 

save(['analyzed_passive_' num2str(sPath(end), '%04i')], 'avgStimResponse', 'avgTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')



%% 
end
