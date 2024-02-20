%% check photobleaching effect within session across different days 
clear; clc; close all
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

fPath = 'Y:\haider\Data\WFI\Mapping\Animals'; 
Animal = 'M221118_1';% use dir
expParadigm = 'PBtest';
expID = {['02-Feb-2023_1'] ['03-Feb-2023']};
sPath = ['Y:\haider\Data\analyzedData\EKK\WFI' filesep Animal filesep expParadigm];
if ~isfolder(sPath)
    mkdir(sPath)
end 
ledIntensity = 10; %in mW 
imDuration = 13; % in min
isi = [20 10 5 1];
% xpix = 150:217; ypix = 172:282;
violet = [0.78,0.65,0.80];
blue = [0.61,0.83,0.98 ];
%% selectinng ROIs
cPath =[fPath filesep Animal filesep expParadigm filesep expID{1,1}];
imshow(imread([cPath filesep 'Snapshot_1.jpg'])); roi = drawrectangle(gca); 
ypix = roi.Position(1):roi.Position(1)+roi.Position(3); 
xpix = roi.Position(2):roi.Position(2)+roi.Position(4);
%%
tic
figure;
for i= 1:length(expID) 
    fclose('all');
    fName = 'Frames'; %name format for imaging files
    dataType = 'uint16'; %type of imaging data usually uint16 or uint8
    cPath =[fPath filesep Animal filesep expParadigm filesep expID{1,i}];
    %% load data
    rawVids = dir([cPath filesep fName '_*']); %video files
    
    % get trial numbers for each file and sort in ascending order
    trials = zeros(1,length(rawVids));
    for x = 1 : length(rawVids)
        [~,a] = fileparts(rawVids(x).name);
        a = textscan(a,'%s','delimiter','_');
        trials(x) = str2double(a{1}{end});
    end
    [trials,sortIdx] = sort(trials,'ascend');
 
    %% photobleaching within experiment
    for j = trials
        cFile = [cPath filesep rawVids(sortIdx(j)).name];
        try
            load([cPath filesep 'frameTimes_' num2str(trials(j), '%04i')], 'imgSize', 'frameTimes') 
            [~, cData] = loadRawData(cFile, 'Frames', [], imgSize);
        catch
            [header, cData] = loadRawData(cFile, 'Frames', dataType); %if file was written was old version, the image size is stored in the binary file
            if isempty(header) %if no header is found frametimes should be
                load([cPath filesep 'frameTimes_' num2str(trials(j))], 'frameTimes');
            else
                frameTimes = header(1:end-4);
            end
        end
        if length(size(cData)) == 4
            if size(cData,3) == 3 % convert rgb image to gray
                cData =  squeeze(sum(cat(3, 0.2989 .* cData(:,:,1,:), 0.5870 .* cData(:,:,2,:), 0.1140 .* cData(:,:,3,:)), 3));
            else
                cData = squeeze(cData(:,:,1,:));
            end
        end
        % dont use first and last frames
        cData = cData(:,:,31:end-30);
        frameTimes = frameTimes(31:end-30);
        
        allData{j} = mean(reshape(cData(floor(xpix),floor(ypix),:),[], size(cData,3),1))'; %keep frame average
        allTimes{j} = frameTimes;
        
        temp = zscore(single(allData{j}));
        viol = find(temp < 0);
        bInd = true(1,length(temp));
        bInd(viol) = false;
%         allData_v{i,j} = allData{j}(~bInd,1);
%         allTimes_v{i,j} = allTimes{j}(~bInd,1);
        allData_b{i,j} = allData{j}(bInd,1);
        allTimes_b{i,j} = allTimes{j}(bInd,1);
        allData_dff{i,j} = normDFF (allData_b{i,j}, 1:length(allData_b{i,j}));
%         plot((allTimes_b{i,j}-allTimes_b{i,j}(1)).*1440, allData_b{i,j}, 'Color', (blue-0.1*j)); hold on
%         plot((allTimes_v{i,j}-allTimes_v{i,j}(1)).*1440, allData_v{i,j},  'Color', (violet-0.1*j)); 
        disp(['done loading session ' num2str(j), '/' num2str(length(trials))])
    end
%     legend(['Blue - Session 1'], ['Violet - Session 1'], ['Blue - Session 1'],...
%         ['Violet - Session 2'],['Blue - Session 2'],...
%         ['Violet - Session 3'],['Blue - Session 3'],...
%         ['Violet - Session 4'],['Blue - Session 4'], ['Blue - Session 5'], ['Violet - Session 5']);
%     xlabel('Time (minutes)'); ylabel('Mean fluorescence'); axis square
    %%
    % combine all trials
    allData = cat(1, allData_b{i,:}); %concatenate all sessions' blue data
    allData_dff1 = cat(1, allData_dff{i,:});
    allTimes = cat(1, allTimes_b{i,:});
    allTimes = (allTimes - allTimes(1)) .* 1440; %convert to minutes
    allData(zscore(diff(allTimes)) > 1) = NaN; %to avoid lines between sessions (end of session and start of ISI)
    allData_dff1(zscore(diff(allTimes)) > 1) = NaN;
    %     allData1(~isnan(allData_b)) = smooth(allData_b(~isnan(allData_b)), 5); % do some smooth (moving average filter) 300ms binning
    allData(~isnan(allData)) = smooth(allData(~isnan(allData)), 50); % do some smooth (moving average filter) 300ms binning
    allData_dff1(~isnan(allData_dff1)) = smooth(allData_dff1(~isnan(allData_dff1)), 50);
    
    s= subplot(2,5,i);
    plot(allTimes,allData, 'linewidth', 1, 'Color', blue-0.3);
    xlabel('Time (minutes)'); ylabel('Mean fluorescence'); axis square
    ylim(s, [1700 2600]); 
    cPath1 = strsplit(cPath,filesep);
    title(sprintf('Mean frame intensity\n %s/%s/%s', cPath1{1,end-2:end}));
    s1 = subplot(2,5,i+5);
    plot(allTimes,allData_dff1, 'linewidth', 1, 'Color', blue-0.3);
    xlabel('Time (minutes)'); ylabel('df/f'); axis square
    ylim(s1, [-0.1 0.2]);
    clear allData allTimes
    disp(['done loading experiment ' expID{i}])

end 
toc 
%%
if ~isempty(sPath)
    savefig([sPath filesep Animal 'PhotobleachingAssessment_days.fig']); %EK
    saveas(gcf, [sPath filesep Animal 'PhotobleachingAssessment_days.png']); %EK
end
%% 
bleaching.max = cellfun(@max, allData_b);
bleaching.min = cellfun(@min, allData_b);
bleaching.mean = cellfun(@mean, allData_b);
bleaching.bleach_rate_session = ((bleaching.max - bleaching.min)./bleaching.max./imDuration) *100; % percentage/min
bleaching.bleach_rate_exp = ((bleaching.max(:,1) - bleaching.min(:,5))./bleaching.max(:,1)) *100; %percentage
[exp,session] = size(allData_b);
%get how much the fluorescence was recovered from the previous session
for i = 1: exp
    for j = 1:session-1      
        bleaching.recover_fluor(i,j)= (mean(allData_b{i,j+1}(1:10))-mean(allData_b{i,j}(end-9:end))); 
    end
end 

save([sPath filesep 'bleaching'], 'bleaching')

figure('Name', 'max fluorescence'); plot(bleaching.max', 'o:', 'linewidth', 2); 
xlabel('session'); ylabel('max fluorescence'); 
Legend = strsplit(sprintf('day%d ',1:5)); 
legend(Legend{1:5})
savefig([sPath filesep 'max fluorescence.fig'])
saveas(gcf, [sPath filesep 'max fluorescence.png']);
figure('Name', 'min fluorescence'); plot(bleaching.min', 'o:','linewidth', 2); 
xlabel('session'); ylabel('min fluorescence'); 
legend(Legend{1:5})
savefig([sPath filesep 'min fluorescence.fig'])
saveas(gcf, [sPath filesep 'min fluorescence.png']);

figure; 
subplot (1,2,1)
plot(isi, bleaching.recover_fluor', 'o:', 'linewidth', 1.5)
set(gca, 'xdir','reverse'); xlim ([0 21]); legend(Legend{1:5});
xlabel('inter-session interval (min)'); ylabel('recovered fluorescence');
axis square
subplot (1,2,2) 
plot(bleaching.bleach_rate_session', 'o:', 'linewidth', 1.5)
xlim([0 6])
xlabel('session'); ylabel('bleaching rate (%/min)');
legend(Legend{1:5}); axis square
savefig([sPath filesep 'bleaching assessment.fig'])
saveas(gcf, [sPath filesep 'bleaching assessment.png']);