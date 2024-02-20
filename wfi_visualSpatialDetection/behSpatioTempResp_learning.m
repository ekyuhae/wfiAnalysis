function [avgStimResponse, avgTempTrace, binocOn, monocOn] = behSpatioTempResp_learning(paramSorted, allData, bFrameTimes, sPath, fwindow,img,type, sessNr, expID_ret, hva, flag)
% get temporal trace and spatial activity map from the imaging during spatial detection
% task experiment.
% modified version of passive grating analysis (barmapFullcont_gratings) 
% EK Feb23
% NA 6/1/2023 - added in hva input and change in plots based on this
%used in [1] wfiBehavior
% EK 8/31 - modified to analyize only one location (before block switch
% phase)
%% 

for i = 1:length(paramSorted)
    contrasts(i) = paramSorted{1,i}(2);
end 
% idx = find(contrasts == 0, 1, 'last');
bcont = contrasts; %was indexing out of bounds?
% mcont = contrasts(idx:end);
binocc = paramSorted;
% binocc = paramSorted(1:idx-1); %was indexing out of bounds?
% monocc = paramSorted(idx:end);

%% WFI barmapping analysis

binocOn = [];
% monocOn = [];
binocFrames = [];
binocTemp = [];
% monocFrames = []; 
% monocTemp = [];
fr = img.sRate;

avgStimResponse.binoc = cell(1,length(bcont));
% avgStimResponse.monoc = cell(1,length(mcont));
indivTempTrace.binoc = cell(1,length(bcont));
% indivTempTrace.monoc = cell(1,length(mcont));
avgTempTrace.binoc = zeros(length(bcont), fr*3+1);
% avgTempTrace.monoc = zeros(length(mcont), fr*3+1);

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

    if paramSorted{1,1}(1) >= 45 
        stimloc = 'monoc';
    else 
        stimloc = 'binoc';
    end 

    % temporal analysis
    if cont == 1 && sessNr ==1 && strcmp(type, 'whole') %
        if ~isempty(hva)
            sPath = fullfile(sPath, '..');
        end
        if exist(fullfile(sPath, 'coordinates.mat'), 'file') == 2
            load(fullfile(sPath, 'coordinates.mat'));
            disp(['coordinates.mat has been loaded']);
        else
            cord = registerRetMap(sPath, 2, expID_ret, [], flag); % for behavior exp (2), register retinotopy
        end
        
        % ROI analysis; select ROI for the first time
        if ~isempty(binocTemp)
            [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData, [], cord, sPath, stimloc);
            avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
        end
    else
        if ~isempty(hva) 
            sPath= fullfile(sPath, '..');
            cfile =[sPath '\..\' num2str(1) filesep hva filesep 'analyzed_active_1_whole' hva '.mat']; % EK added 'whole'
        else
            cfile =[sPath(1:end-1) num2str(1) filesep 'analyzed_active_1_whole.mat'];
        end
        if exist(cfile) % for second session and so on, automatically load roi that was chosen before and get avg roi trace (only if stim location is same as previous one)
            load(cfile,'roi', 'cord')
            prevLoc =  load(cfile, 'paramSorted');
            % disp('ROI has been loaded')
            if ~isempty(binocTemp)
                if (~isfield(roi, stimloc) || ~isequal(prevLoc.paramSorted{1,1}(1), paramSorted{1,1}(1))) && cont ==1 && strcmp(type,'whole') %EK OCT23
                    [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,[], cord, sPath, stimloc);
                else
                    [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,roi.binoc, cord, sPath,stimloc);
                end
                avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
            end
        else
            if ~isempty(binocTemp)
                if ~isfield(roi, stimloc)
                    [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,[], cord, sPath, stimloc);
                else
                    [indivTempTrace.binoc{cont}, roi.binoc] = getIndividTrace(binocTemp,btraceLength, allData,roi.binoc, cord, sPath, stimloc);
                end
                avgTempTrace.binoc(cont,:) = mean(indivTempTrace.binoc{cont},1);
            end
        end
    end
  
        binocFrames = [];
        binocTemp = [];
        if ~isempty(hva)
            sPath = sPath(1:end-2);
        end
end



cd(sPath)
% plotAvgMap_beh(contrasts, avgStimResponse,cord, sPath, '', hva) % contrast range for binoc and monoc is different 5/11 fix ity

sessDir = fileparts(pwd); 
% if ~strcmp(sessDir, [sessDir filesep num2str(1)]) && length(sessDir) > length([sessDir(1:end-2) '\' num2str(1)]) %if it's subfolder
%   sessDir = [sessDir '\..'];
% end

% if ~isempty(hva) || reactionTriggered
%     sessDir = [sessDir '\..'];
% end
    load([sessDir filesep num2str(1) filesep 'contours.mat']) %
    load([sessDir filesep num2str(1) filesep 'coordinates.mat'], 'snap')
contours = imread([sessDir filesep num2str(1) filesep 'overlaid_contourMap.png']);
% catch
%     load([sPath filesep 'contours.mat'])
%     load([sPath filesep 'coordinates.mat'])
%     contours = imread([sPath filesep 'overlaid_contourMap.png']);
% % end
contours = rgb2ind(contours, 256);

%%
% idx = find(contrasts == 0, 1, 'last');

% bcont = contrasts(1:idx);
for i = 1: length(avgStimResponse.binoc)
    if ~isempty(avgStimResponse.binoc{i})
        binoc{i} = avgStimResponse.binoc{i};
    else
        binoc{i}= NaN;
    end
end

m = figure;
% colorscale = [min(min(min(binoc))) max(max(max(binoc)))];
% tiledlayout(1, length(contrasts),'TileSpacing', 'tight', 'Padding', 'compact'); % show avg stim response activity for each white bar location
tiledlayout(1, length(bcont),'TileSpacing', 'compact', 'Padding', 'compact'); % show avg stim response activity for each white bar location
for i = 1: length(bcont)
    binoc_resized{i} = imresize(binoc{i}, [size(contours,1), size(contours,2)]);
    [binoc_resized{i},~] = imalign_s(binoc_resized{i}, contours, cord.moved,cord.fixed);
    
    nexttile
    imagesc(binoc_resized{i}); colormap viridis; freezeColors; colorbar
    hold on
    contour(azi, azi_levels, 'k', 'LineWidth', 0.02); daspect([1 1 1]);
    contour(nanmean(imout,3), 'k', 'LineWidth',0.02); axis ij; axis off
    axis image
    
    title(sprintf('%d%%', bcont(i)*100));
end
if isempty(hva)
    sgtitle (['avg activity map - binoc ' type ' trials']);
else
    sgtitle(['avg activity map - ' hva(1:2) ' ' type ' trials'])
end

set(m, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 40 5])
if isempty(hva)
    savefig(m,  [sPath filesep 'avgActivityMap_binoc_' type 'Trials.fig']);
else 
    savefig(m, [sPath filesep 'avgActivityMap_' hva(1:2) '_' type 'Trials.fig'])
end
% saveas(m,  ['avgActivityMap_binoc_' type 'Trials.png']);
% if isempty (hva)
% plotAvgMap_beh(contrasts, avgStimResponse,cord, sPath, [], hva) % contrast range for binoc and monoc is different 5/11 fix ity
% else
% 
%    plotAvgMap_beh(contrasts, avgStimResponse,cord, fullfile(sPath, '..'), [], hva) % contrast range for binoc and monoc is different 5/11 fix ity 
% end 
k = linspace(0.8, 0, length(bcont));
fr = fr/2;
t = 0:1/fr:(length(avgTempTrace.binoc)-1)/fr; % convert frame into timescale 
t = t - t(fr*2+1); 

b= figure;
subplot (1,2,1)
legend_names = cell(1, length(bcont));
for i = 1:length(bcont)
    legend_names{i} = sprintf('%d%%', bcont(i)*100);
    color = [k(i), k(i), k(i)];    
%     y_adj = imadjust(avgTempTrace.binoc(i,:), [], [], i/length(contrasts));
%     plot(t, y_adj, 'Color', color);
     plot(t, avgTempTrace.binoc(i,:), 'Color', color);
    hold on;

end
xline(t(fr*2+1), ':k') ;
if isempty(hva)
    title(['binoc ROI trace - ' type ' trials']);
else
    title([hva(1:2) ' ROI trace - ' type ' trials'])
end
xlabel('time (s)'); ylabel('dF/F_0')
legend(legend_names); legend('boxoff'); axis square
% if avgTempTrace.binoc ~= 0
%     yLim = [min(min(avgTempTrace.binoc)) max(max(avgTempTrace.binoc))];
%     ax = findobj(b,'type','axes'); % Find all axes handles in the figure
%     set(ax,'ylim',yLim);
% end
hold off

if ~isempty(type)
    savefig(b, ['ROI traces_' type ' trials.fig'])
%     saveas(b, ['ROI traces_' type ' trials.png'])
%     savefig(m, ['monoc ROI trace_' type ' trials.fig'])
%     saveas(m, ['monoc ROI trace_' type ' trials.png'])
else
    savefig(b, 'ROI traces_whole trials.fig')
%     saveas(b, 'ROI traces_whole trials.png')
%     savefig(m, 'monoc ROI trace_whole trials.fig')
%     saveas(m, 'monoc ROI trace_whole trials.png')

end

if isempty(hva)
% save(['analyzed_active_' num2str(sessNr) '.mat'], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')

    save(['analyzed_active_' num2str(sPath(end), '%04i') '_' type], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
else
    save(['analyzed_active_' num2str(sessDir(end)) '_' hva '_' type], 'avgStimResponse', 'avgTempTrace', 'indivTempTrace', 'roi', 'contrasts', 'paramSorted', 'cord')
end

%% 
end