% EK May23
clc; clear; close all

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

Animal = 'M230220_1';
expID = '22-May-2023'; % this is WFI experiment ID format
expID_ret = '04-Apr-2023_1'; % for retinotopy registration and identifying the visual areas 
% expID_ret = '16-Mar-2023';
singleSession = false;
sessNr = []; % enter the imaging session number(X) from Frames_2_xxx_xxx_uint16_000X 
cPath = [Animal filesep 'behavior' filesep expID filesep num2str(sessNr)];
savePath = (['Y:\haider\Data\analyzedData\EKK\WFI\' cPath]); % enter path to saving the analyzed data
type = 'hit'; %  behavior performance type
% expInfoNr = '1'; % enter the corresponding expInfo session number (does not always match with WFI sessNr) 
date = '2023-05-22'; 
[~, totalTrials, ~, beh_all, activeOut] = attention_effects_EK(Animal, date, {'1','2','3','4'});

%% make movie of response in single session
if singleSession
    cd(savePath)
    cFile = "analyzed_active_1_hit.mat";
    if exist(cFile)
        load(cFile, 'paramSorted')
    end

    for i = 1:length(paramSorted)
        contrasts(i) = paramSorted{1,i}(2);
    end
    idx = find(contrasts == 0, 1, 'last');
    bcont = contrasts(1:idx-1);
    mcont = contrasts(idx:end);
    binoc = paramSorted(1:idx-1);
    monoc = paramSorted(idx:end);

    temp = [];
    binocOn = [];
    monocOn = [];
    binocTemp = [];
    monocTemp = [];
    datab =[]; datam =[];
    fr = img.sRate;

    for cont = 1:length(bcont)
        binocOn{cont} = findStimOnFrame(binoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(binocOn{cont})
            for i = 1:length(binocOn{cont})
                binocTemp= cat(2,binocTemp, binocOn{cont}(i)-fr:binocOn{cont}(i)+fr*2); % 6s window - 2s prestim, 4s poststim
            end
            btraceLength = length(binocTemp) / length(binocOn{cont});
            binocTemp = reshape(binocTemp,btraceLength,[])';
            for t= 1: size(binocTemp,1)
                for frame = 1:size(binocTemp,2)
                    temp(:,:,frame,t) = allData(:,:, binocTemp(t,frame));
                end
            end
            dSize = size(temp);
            datab = cat(4, datab, temp);
            dSize_b = size(datab);
        end
        binocTemp = [];
        temp = [];
    end
    bAvgMovie = mean(datab,4);
    compareMovie(bAvgMovie);

    for  cont = 1:length(mcont)
        monocOn{cont} = findStimOnFrame(monoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(monocOn{cont})
            for i = 1:length(monocOn{cont})
                monocTemp= cat(2,monocTemp, monocOn{cont}(i)-fr:monocOn{cont}(i)+fr*2);
            end
            mtraceLength = length(monocTemp) / length(monocOn{cont});
            monocTemp = reshape(monocTemp,mtraceLength,[])';
            for t= 1: size(monocTemp,1)
                for frame = 1:size(monocTemp,2)
                    temp(:,:,frame,t) = allData(:,:, monocTemp(t,frame));
                end
            end
            dSize = size(temp);
            datam = cat(4, datam, temp);
            dSize_m = size(datam);
        end
        monocTemp = [];
        temp = [];
    end
    mAvgMovie = mean(datam,4);
    compareMovie(mAvgMovie);
    clearvars -except allData img stimOn bFrameTimes vFrameTimes
end

%% make movie of response in multiple sessions
cPath = [Animal filesep 'behavior' filesep expID];
savePath = (['Y:\haider\Data\analyzedData\EKK\WFI\' cPath]); % enter path to saving the analyzed data

cd(savePath)
sessions = dir(savePath);
sessions = sessions([sessions(:).isdir]); % keep only directories
sessions = sessions(3:end);
count = 0;
for i = 1:length(sessions)
    if regexp(sessions(i).name, '^\d+$')
        count = count + 1;
    end
end
nrSessions= count;

%%
datab =[]; datam =[];
for sess = 1: nrSessions
    load([savePath filesep num2str(sess) filesep 'preprocessed_WFIdata_' num2str(sess, '%04i') '.mat'], 'allData', 'bFrameTimes', 'img', 'stimOn')
    cFile = [savePath filesep num2str(sess) filesep 'analyzed_active_' num2str(sess) '_' type '.mat'];
    if exist(cFile)
        load(cFile, 'paramSorted')
    end

    for i = 1:length(paramSorted)
        contrasts(i) = paramSorted{1,i}(2);
    end
    idx = find(contrasts == 0, 1, 'last');
    bcont = contrasts(1:idx-1);
    mcont = contrasts(idx:end);

    binoc = paramSorted(1:idx-1);
    monoc = paramSorted(idx:end);

    temp = [];
    binocOn = [];
    monocOn = [];
    binocTemp = [];
    monocTemp = [];
    
    fr = img.sRate;

    for cont = 1:length(bcont)
        binocOn{cont} = findStimOnFrame(binoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(binocOn{cont})
            for i = 1:length(binocOn{cont})
                binocTemp= cat(2,binocTemp, binocOn{cont}(i)-fr:binocOn{cont}(i)+fr*2); % 6s window - 2s prestim, 4s poststim
            end
            btraceLength = length(binocTemp) / length(binocOn{cont});
            binocTemp = reshape(binocTemp,btraceLength,[])';
            for t= 1: size(binocTemp,1) % timepoint for corresponding stim on frame
                for frame = 1:size(binocTemp,2)
                    temp{cont}(:,:,frame,t) = allData(:,:, binocTemp(t,frame));
                end
            end
            dSize = size(temp{cont});
            datab = cat(4, datab, temp{cont});
            dSize_b = size(datab);
        end
        binocTemp = [];
    end
    temp = [];
    binoc = [];

    for  cont = 1:length(mcont)

        monocOn{cont} = findStimOnFrame(monoc{cont}(3:end), 50, bFrameTimes);
        if ~isempty(monocOn{cont})
            for i = 1:length(monocOn{cont})
                monocTemp= cat(2,monocTemp, monocOn{cont}(i)-fr:monocOn{cont}(i)+fr*2);
            end
            mtraceLength = length(monocTemp) / length(monocOn{cont});
            monocTemp = reshape(monocTemp,mtraceLength,[])';
            for t= 1: size(monocTemp,1)
                for frame = 1:size(monocTemp,2)
                    temp{sess}(:,:,frame,t) = allData(:,:, monocTemp(t,frame));
                end
            end
            dSize = size(temp);
            datam = cat(4, datam, temp{sess});
            dSize_m = size(datam);

        end
        monocTemp = [];

    end
        temp = [];
        monoc =[];
end
    bAvgMovie = mean(datab,4);
    mAvgMovie = mean(datam,4);
%%
load([savePath filesep num2str(1) filesep 'contours.mat'])
load([savePath filesep num2str(1) filesep 'coordinates.mat'])
contours = imread([savePath filesep num2str(1) filesep 'overlaid_contourMap.png']);
contours = rgb2ind(contours, 256);

for i = 1:size(bAvgMovie,3)
    bAvgMovie_resized(:,:,i) = imresize(bAvgMovie(:,:,i), [size(contours,1), size(contours,2)]);
    [bAvgMovie_resized(:,:,i),~] = imalign_s(bAvgMovie_resized(:,:,i), contours, cord.moved,cord.fixed);
end 

for i = 1:size(mAvgMovie,3)
    mAvgMovie_resized(:,:,i) = imresize(mAvgMovie(:,:,i), [size(contours,1), size(contours,2)]);
    [mAvgMovie_resized(:,:,i),~] = imalign_s(mAvgMovie_resized(:,:,i), contours, cord.moved,cord.fixed);
end 
save(['3dStack_' type '.mat'], "bAvgMovie_resized","mAvgMovie_resized","datam","datab");
compareMovie(bAvgMovie_resized);
compareMovie(mAvgMovie_resized);