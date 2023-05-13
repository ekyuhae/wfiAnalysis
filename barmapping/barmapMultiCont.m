function [avgStimResponse, avgTempTrace] = barmapMultiCont(paramSorted, allData, bFrameTimes, sPath, window, img, multicont, sessNr, expID_ret)
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
for cont = 1:length(cont1) % get response maps from different contrasts (contrast sorted in descending order)
    if cont1(cont) ~= 0
        for locc = 1:length(loc) % locations
            %% spatial analysis
            frameSpatio = [];
            frameTemp =[];
            
            barFrames = findStimOnFrame(contrasts{cont}(:,locc),50,bFrameTimes);             
            for i = 1:length(barFrames) % should be same as number of bar reapeats
                if barFrames(i) ~= 0
                frameSpatio = cat(2,frameSpatio, barFrames(i):barFrames(i)+window); % get response window frames from sitm onset and take average
                frameTemp= cat(2,frameTemp, barFrames(i)-fr:barFrames(i)+fr); % for the temporal analysis, save frames for +/- window from stimulus onset
                end
            end 
            avgStimResponse{cont}(:,:,locc) = nanmean(allData(:, :, frameSpatio), 3);
            
           if cont1(cont) ~= 0 && cont == 1 && sessNr ==1 % 1st contrast, 1st session, get all rois according to different locations         
               if locc==1
                   if exist(fullfile(sPath, 'coordinates.mat'), 'file') == 2
                       load(fullfile(sPath, 'coordinates.mat'));
                       disp(['coordinates.mat has been loaded']);
                   else
                       cord = registerRetMap(sPath, 1, expID_ret);
                   end
               end 
               [avgTempTrace{cont,locc}, roi{locc}] = getAvgTemporalResponse(frameTemp,img, allData,[], cord, sPath); %binoc location; 0.6 deg              
                
           else
               cfile =[sPath(1:end-1) num2str(1) filesep 'analyzed_barmapping_1.mat'];
% cfile = 'Y:\haider\Data\analyzedData\EKK\WFI\M221118_1\barmapping\09-Mar-2023\1\3ROIs\analyzed_barmapping_s.mat';
               if exist(cfile)
                   load(cfile,'roi', 'cord')
               end 
               [avgTempTrace{cont,locc},~] = getAvgTemporalResponse(frameTemp,img, allData,roi{locc}, cord, sPath);
           end 
            
        end        
    else % for zero contrast (control trial)
         frameSpatio = [];
         frameTemp =[];
         barFrames = findStimOnFrame(contrasts{cont},50,bFrameTimes);
         for i = 1:length(barFrames) % should be same as number of bar reapeats
             frameSpatio = cat(2,frameSpatio, barFrames(i):barFrames(i)+window); % get response window frames from sitm onset and take average
             frameTemp= cat(2,frameTemp, barFrames(i)-fr:barFrames(i)+fr); % for the temporal analysis, save frames for +/- window from stimulus onset
         end
         avgStimResponse{cont} = nanmean(allData(:, :, frameSpatio), 3);
         [avgTempTrace{cont,1},~] = getAvgTemporalResponse(frameTemp,img, allData,roi{1}, cord, sPath);
   
    end
end 

%%
location.loc = loc; 
location.loc_gr = loc_gr;
% contrasts = cont1;
for cont = 1:length(cont1)
    if cont1(cont) ~= 0
        plotAvgMap_bar(cont1(cont), avgStimResponse{cont}, loc, cord, sPath);
    else
        plotAvgMap_bar(cont1(cont), avgStimResponse{cont}, loc_gr, cord, sPath);
    end
end
plotTempTrace_bar(avgTempTrace, cont1, img.sRate, location, multicont);

save(['analyzed_barmapping_' num2str(sPath(end), '%04i')], 'avgStimResponse', 'avgTempTrace', 'roi', 'cont1','location', 'cord')
