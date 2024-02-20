function [avgStimResponse, avgTempTrace] = barmapMultiCont_lockROI(paramSorted, allData, bFrameTimes, sPath, window, img, multicont, sessNr, expID_ret, flag)
% bar mapping full contrast analysis. updated version of WFIbarmappingAnalysis2
% plot activity maps corresponding to different bar locations
% EK Feb23
%% 
% for multi contrast - combine b&w from high contrast to low contrast 
color = unique(paramSorted(1,:));
for i = 1: floor(length(color)/2)+1
    b = paramSorted(3:end, paramSorted(1,:)==color(i));
    w = paramSorted(3:end, paramSorted(1,:) == color(end-(i-1)));
    if ~isequal(b,w)
        contrasts{i} = vertcat(b,w);
    else 
        contrasts{i} = b;
    end
end 
%% WFI barmapping analysis 
for i = 1:length(contrasts)
    avgStimResponse{i}= zeros(size(allData,1), size(allData,2), size(contrasts{i},2));
end
cont1 =abs(100*(color - color(color==128))/abs(color(color==128))); % convert color into relative contrast change
cont1 = cont1(1:find(color==128));
% loc = paramSorted(2,find(paramSorted(1,:)==0));
loc = paramSorted(2,1:17); % change accordingly if the number of location is not 17 
loc_gr = paramSorted(2,find(paramSorted(1,:)==128));

fr = img.sRate;
% fr = floor(img.sRate*0.25); % for 30Hz, give +- fr = img.sRate frames, 
cfile =[sPath(1:end-1) num2str(1) filesep 'analyzed_barmapping_1.mat'];
if exist(cfile)
    load(cfile,'roi', 'cord')
end

for cont = 1:length(cont1) % get response maps from different contrasts (contrast sorted in descending order)
    if cont1(cont) ~= 0
        for ii = 1: length(roi)
            for locc = 1:length(loc) % locations
                %% spatial analysis
                frameSpatio = [];
                frameTemp =[];

                barFrames{locc} = findStimOnFrame(contrasts{cont}(:,locc),50,bFrameTimes);
                for i = 1:length(barFrames{locc}) % should be same as number of bar reapeats
                    if barFrames{locc}(i) ~= 0
                        % frameSpatio = cat(2,frameSpatio, barFrames(i):barFrames(i)+window); % get response window frames from sitm onset and take average
                        frameTemp= cat(2,frameTemp, barFrames{locc}(i)-fr:barFrames{locc}(i)+fr); % for the temporal analysis, save frames for +/- window from stimulus onset
                    end
                end

                traceLength = length(frameTemp)/length(barFrames{locc});

                % if cont1(cont) ~= 0 && cont == 1 && sessNr ==1 % 1st contrast, 1st session, get all rois according to different locations
                if locc==1 && ii ==1
                    if exist(fullfile(sPath, 'coordinates.mat'), 'file') == 2
                        load(fullfile(sPath, 'coordinates.mat'));
                        disp(['coordinates.mat has been loaded']);
                    end
                end

                [avgTempTrace{cont,ii}(:,locc),~] = getAvgTemporalResponse(frameTemp,traceLength, allData,roi{ii}, cord, sPath, []); %binoc location; 0.6 deg
                disp(['roi ' num2str(ii) ' done']);
            end
        end
    end
end

for i = 1: length(roi) % ROIs
    for ii = 1:length(loc) % bar location
        tmp(i,ii) = getPeakPostStim(avgTempTrace{i}(:,ii)', 15, 0.5,2,0);
        % tmp(i,ii) = mean(avgTempTrace{i}(31:36,ii),1);
    end 
end 

figure;

for i = 4:length(loc)
plot(loc, tmp(i,:)); hold on
end 

for i = 1:length(roi)
    [~,idx(i)] = max(tmp(i,:));
end 
roi_barloc = loc(idx);
%%
% location.loc = loc; 
% location.loc_gr = loc_gr;
% % contrasts = cont1;
% % for cont = 1:length(cont1)
% %     if cont1(cont) ~= 0
% %         plotAvgMap_bar(cont1(cont), avgStimResponse{cont}, loc, cord, sPath);
% %     else
% %         plotAvgMap_bar(cont1(cont), avgStimResponse{cont}, loc_gr, cord, sPath);
% %     end
% % end
% plotTempTrace_bar(avgTempTrace, cont1, img.sRate, location, multicont);
cd(sPath)
save(['tuningCurveResp.mat'], 'avgTempTrace', 'tmp', 'loc', 'roi_barloc')
