 %% WFI Behavior Experiment Pipeline - Analysis for ATTENTIONAL EFFECTS
% analysis after running wfiBehavior (analysis for visual detection; df/f
% PSTHs for different tial types) since it requires
% 'analyed_active_(sessNr).mat' for pulling ROI
% pull the data from the session that has attentional signs in binoc and/or
% monoc blocks
% 
% EK Oct 23
% EK 12/8 edited
% EK Feb 24 edited; HVA option added

clc; clear; close all

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

Animal = 'M230831_1';
% expID = {'17-Aug-2023' '21-Aug-2023' '23-Aug-2023' '01-Sep-2023' '04-Sep-2023' '26-Sep-2023'...
%     '27-Sep-2023' '03-Oct-2023' '04-Oct-2023' '17-Oct-2023', '24-Oct-2023'}; % this is WFI experiment ID format

hva = ''; %higher visual areas being analyzed (i.e. 'RL_AM' or 'LM_PM' or '' for V1)
a.binoc ={[1]; [1 3]; [4]; []; [3]; [6 7]; []; [5 6]; [2 3 4]; [1 2 4]; []};
a.monoc = {[]; []; []; [2]; []; []; [9 10 11]; []; []; [1 2 3]; [1 2 3]} ; % WFI session Number 

b.binoc = {[1]; [1 3]; [4]; []; [3]; [6 7]; []; [5 6]; [2 3 4]; [1 2 4]; []};
b.monoc = {[]; []; []; [2]; []; []; [9 10 11]; []; []; [1 2 3]; [5 6 7]} ; % expInfo session Number
sessNr = [a.binoc, a.monoc];
expInfoNr = [b.binoc, b.monoc];
clear a b

cd(['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior'])
if ~isfolder('attention_over_days')
    if isempty(hva)
        mkdir('attention_over_days')
    else 
        mkdir(['attention_over_days_' hva])
    end 
end

cd('attention_over_days')
if ~isfile ('key_attn.mat')
    key.Animal = Animal;
    key.expID = expID;
    key.expInfoNr = expInfoNr;
    key.sessNr = sessNr;
    save('key_attn', "key");
else 
    % load("key_attn.mat");
    load("key_attn_updated.mat");
    expID = key.expID;
    expInfoNr = key.expInfoNr;
    sessNr = key.sessNr; 
end

%% pulling the analyzed mat file (output from running wfiBehavior) 
for cDay = 1:length(expID) % CHANGE TO 1
    date = datestr(expID{cDay}(1:11), 'yyyy-mm-dd'); % expInfo date
    cSess = unique(cat(2,sessNr{cDay,1}, sessNr{cDay,2}));
    cExpInfoNr = unique(cat(2,expInfoNr{cDay,1}, expInfoNr{cDay,2}));
    binocTemp = [];  monocTemp = [];
    for s = 1: length(cSess)
        expPath = ['Y:\haider\Data\Behavior\expInfo\' Animal filesep date filesep num2str(cExpInfoNr(s))];
        sessID = [date '_' num2str(cExpInfoNr(s)) '_' Animal];
        load ([expPath filesep sessID '_Block.mat']);    
        [out, beh_all, activeOut, trTable] = attention_effects_hitAndTrialAligned(Animal, date, {num2str(cExpInfoNr(s))}); % 8/29 EK changed to attention code from regular analysis code
        close all
        fPath ='Y:\haider\Data\analyzedData\EKK\WFI\';
        cd([fPath filesep Animal filesep 'behavior' filesep expID{cDay} filesep num2str(cSess(s))])
        if s ==1
            if isempty(hva)
                load(['analyzed_active_' num2str(cSess(s)) '.mat'], 'roi', 'cord');
            else
                cd(hva)
                load(['analyzed_active_' num2str(cSess(s)) '_' hva '.mat'], 'roi', 'cord');
                cd(fileparts(pwd))
            end
        end
        load(['preprocessed_WFIdata_' num2str(cSess(s), '%04i') '.mat'], 'allData', 'bFrameTimes', 'vFrameTimes', 'stimOn', 'img');
        fr = img.sRate;
        idx = [ismember(cSess(s), sessNr{cDay,1}) ismember(cSess(s), sessNr{cDay,2})]; % indexing binoc or monoc attention

            %% find the frametime that corresponds to first 1~10 th trials in the block that had attentional signs for each session
        if idx == [1 1] % the session had attn signs for both blocks

            % sort hit trials over each block (e.g. 1,2,3.. 10th hit trials across 3 binoc blocks) and align their stimulus timepoints to frametime
            [bHitTimes, mHitTimes] = behSortParameter_attn(trTable, block, bFrameTimes, vFrameTimes, stimOn);
            bTrNr(s) = max(trTable.trNumInBlock (trTable.boolBinoc));
            mTrNr(s) = max(trTable.trNumInBlock (trTable.boolMonoc));

            % find frame time that corresponds to each stim onset timepoint
            indivTempTrace(s).binoc = cell(1,length(bHitTimes));
            indivTempTrace(s).monoc = cell(1,length(mHitTimes));
            avgTempTrace(s).binoc = zeros(length(bHitTimes), fr*3+1);
            avgTempTrace(s).monoc = zeros(length(mHitTimes), fr*3+1);

            binocOn = cellfun(@(x) findStimOnFrame(x, 50, bFrameTimes), bHitTimes, 'UniformOutput', false);
            monocOn = cellfun(@(x) findStimOnFrame(x, 50, bFrameTimes), mHitTimes, 'UniformOutput', false);

            tmp = cellfun(@(x) arrayfun(@(x) x - fr : x +fr*2, x, 'UniformOutput', false)', binocOn, 'UniformOutput', false);
            binocTemp = cellfun(@(x) cell2mat(x), tmp, 'UniformOutput', false);

            tmp1 = cellfun(@(x) arrayfun(@(x) x - fr : x +fr*2, x, 'UniformOutput', false)', monocOn, 'UniformOutput', false);
            monocTemp = cellfun(@(x) cell2mat(x), tmp1, 'UniformOutput', false);
            clear tmp tmp1 binocOn monocOn

            btraceLength = length(binocTemp{1,1});
            mtraceLength = length(monocTemp{1,1});

            %temporal analysis - requires alignment coordinates and ROI pixels
            sPath= pwd;
            for i = 1: length(binocTemp)
                indivTempTrace(s).binoc{i} = getIndividTrace(binocTemp{i},btraceLength, allData, roi.binoc, cord, sPath, flag);
                avgTempTrace(s).binoc(i,:) = mean(indivTempTrace(s).binoc{i},1);
            end
            for i = 1: length(monocTemp)
                indivTempTrace(s).monoc{i} = getIndividTrace(monocTemp{i},mtraceLength, allData, roi.monoc, cord, sPath, flag);
                avgTempTrace(s).monoc(i,:) = mean(indivTempTrace(s).monoc{i},1);
            end
            binocTemp =[]; monocTemp= []; bHitTimes =[]; mHitTimes =[];

        elseif idx == [1 0] % binoc 
            [bHitTimes, ~] = behSortParameter_attn(trTable, block, bFrameTimes, vFrameTimes, stimOn);
            bTrNr(s) = max(trTable.trNumInBlock (trTable.boolBinoc));

            % find frame time that corresponds to each stim onset timepoint
            indivTempTrace(s).binoc = cell(1,length(bHitTimes));
            avgTempTrace(s).binoc = zeros(length(bHitTimes), fr*3+1);

            binocOn = cellfun(@(x) findStimOnFrame(x, 50, bFrameTimes), bHitTimes, 'UniformOutput', false);

            tmp = cellfun(@(x) arrayfun(@(x) x - fr : x +fr*2, x, 'UniformOutput', false)', binocOn, 'UniformOutput', false);
            binocTemp = cellfun(@(x) cell2mat(x), tmp, 'UniformOutput', false);
            clear tmp binocOn 

            btraceLength = length(binocTemp{1,1}) ;

            sPath= pwd;
            for i = 1: length(binocTemp)
                indivTempTrace(s).binoc{i} = getIndividTrace(binocTemp{i},btraceLength, allData, roi.binoc, cord, sPath, flag);
                avgTempTrace(s).binoc(i,:) = mean(indivTempTrace(s).binoc{i},1);
            end
            binocTemp =[]; bHitTimes =[]; 

        elseif idx == [0 1] % monoc
            [~, mHitTimes] = behSortParameter_attn(trTable, block, bFrameTimes, vFrameTimes, stimOn);
            mTrNr(s) = max(trTable.trNumInBlock (trTable.boolMonoc));

            % indivTempTrace(s).monoc = cell(1,length(mHitTimes));
            % avgTempTrace(s).monoc = zeros(length(mHitTimes), fr*3+1);
            % 
            monocOn = cellfun(@(x) findStimOnFrame(x, 50, bFrameTimes), mHitTimes, 'UniformOutput', false);

            tmp1 = cellfun(@(x) arrayfun(@(x) x - fr : x +fr*2, x, 'UniformOutput', false)', monocOn, 'UniformOutput', false);
            monocTemp = cellfun(@(x) cell2mat(x), tmp1, 'UniformOutput', false);
            clear tmp1 monocOn

            mtraceLength= length(monocTemp{1,1});
            sPath= pwd;
            for i = 1: length(monocTemp)
                indivTempTrace(s).monoc{i} = getIndividTrace(monocTemp{i},mtraceLength, allData, roi.monoc, cord, sPath, flag);
                avgTempTrace(s).monoc(i,:) = mean(indivTempTrace(s).monoc{i},1);
            end
            monocTemp= []; mHitTimes =[];
        end
    end

    %% avg across sessions
    % k = linspace(0.8, 0, length(bcont));
    binocTemp = []; monocTemp =[];

    if isfield(avgTempTrace, 'binoc')
        nrAlignedHits = max(arrayfun(@(x) size(x.binoc,1), avgTempTrace));
        for trial= 1:nrAlignedHits
            for i = 1:length(avgTempTrace)
                if trial <= height(avgTempTrace(i).binoc)
                    binocTemp = [binocTemp; avgTempTrace(i).binoc(trial,:)];
                else
                    binocTemp = binocTemp;
                end
            end
            avgSessionB(trial,:) = mean(binocTemp,1);
            binocTemp =[];
        end
    end
    if isfield(avgTempTrace, 'monoc')
        nrAlignedHits = max(arrayfun(@(x) size(x.monoc,1), avgTempTrace));
        for trial= 1:nrAlignedHits
            for i = 1:length(avgTempTrace)
                if  trial <= height(avgTempTrace(i).monoc)
                    monocTemp = [monocTemp; avgTempTrace(i).monoc(trial,:)];
                else
                    monocTemp = monocTemp;
                end
            end
            avgSessionM(trial,:) = mean(monocTemp,1);
            monocTemp =[];
        end
    end
    cd(fileparts(pwd));
    if ~isfolder('attention_over_sessions')
           mkdir('attention_over_sessions')
    end
    if isempty(hva)
        cd('attention_over_sessions\')
    else
        if ~isfolder(['attention_over_sessions_' hva])
            mkdir(['attention_over_sessions_' hva])
        end 
        cd(['attention_over_sessions_' hva])
    end

    bSess  = sessNr{cDay,1}; mSess = sessNr{cDay,2};
    if exist('mTrNr', 'var') && exist('bTrNr', 'var')
        save(['hit_aligned_dff.mat'],'avgTempTrace', 'bTrNr', 'mTrNr', 'bSess', 'mSess');
    elseif exist('mTrNr', 'var') && ~exist('bTrNr', 'var')
        save(['hit_aligned_dff.mat'],'avgTempTrace', 'mTrNr', 'bSess', 'mSess');
    elseif ~exist('mTrNr', 'var') && exist('bTrNr', 'var') 
        save(['hit_aligned_dff.mat'],'avgTempTrace', 'bTrNr', 'bSess', 'mSess');
    end

    %% plots
    fr = fr/2;
    t = 0:1/fr:90/fr; % convert frame into timescale
    t = t - t(fr*2+1);

    if exist("avgSessionB", "var")
        figure ;
        if height(avgSessionB) < 10
            k = height(avgSessionB);
        else
            k = 10;
        end
        for i = 1: k
            subplot(1,k,i)
            plot(t, avgSessionB(i,:), 'LineWidth', 1.5, 'color','k');
            xline(t(fr*2+1), ':k');
            title(['hit no.' num2str(i)]);
            xlabel('time (s)'); ylabel('dF/F_0')
            axis square; box off
            ylim([-0.05 0.08])
        end
        sgtitle(['Binoc Attentional block: sess ' strsplit(num2str(sessNr{cDay,1}), ' ')])

        savefig(gcf, 'dffAcrossBinocHitTrials')
        saveas(gcf, 'dffAcrossBinocHitTrials.png')
    end

    if exist("avgSessionM", "var")
        figure;
        if height(avgSessionM) < 10
            k = height(avgSessionM);
        else
            k = 10;
        end
        for i = 1: k
            subplot(1,k,i)
            plot(t, avgSessionM(i,:), 'LineWidth', 1.5, 'color','k');
            xline(t(fr*2+1), ':k');
            title(['monoc hit no.' num2str(i)]);
            xlabel('time (s)'); ylabel('dF/F_0')
            axis square; box off
            ylim([-0.05 0.08])
        end
        sgtitle(['Monoc Attentional block: sess ' strsplit(num2str(sessNr{cDay,2}), ' ')])
        savefig(gcf, 'dffAcrossMonocHitTrials')
        saveas(gcf, 'dffAcrossMonocHitTrials.png')
    end
    %%
    if exist("bTrNr", "var")
        bTrNr = max(bTrNr); bPeak = nan(1,bTrNr); bBaseline = nan(1,bTrNr); bBaselineStd = nan(1,bTrNr);
        for i = 1: size(avgSessionB,1)
            [bPeak(i), peak_index] = getPeakPostStim(avgSessionB(i,:), 15, 1, 2, 1);
            bBaseline(i) = mean(avgSessionB(i,1:30),2); 
            % bBaselineStd(i) = std(avgSessionB(i,1:30),0,2);
        end

        figure;
        
        p1 = plot(1:bTrNr, bPeak, 'ko-', 'lineWidth', 1.5); hold on
        p2 = plot(1:bTrNr, bBaseline, 'k:', 'lineWidth', 1.5);
        [r,p,xfit,yfit] = fit_line(1:5, bPeak(1:5)');
        p3 = plot(xfit, yfit, 'r', 'LineWidth', 2);
        legend ([p1 p2], {'postStim', 'baseline'})

        % errorbar (1:bTrNr, bBaseline, bBaselineStd ,'-o', 'CapSize',0, 'color', [0.5 0.5 0.5])
        xlabel('trial number'); ylabel('dF/F_0')
        title('binoc'); axis square; box off;
        legend boxoff
        xlim([0 bTrNr]); ylim([-0.02 0.08]);
        set(findobj(gcf, 'type', 'axes'), 'FontSize', 12)
        savefig(gcf, 'peak_vs_BinocHitTrialNumber')
        saveas(gcf, 'peak_vs_BinocHitTrialNumber.png')
    end

    if exist("mTrNr", "var")
        mTrNr = max(mTrNr); mPeak = nan(1,mTrNr); mBaseline = nan(1,mTrNr); mBaselineStd = nan(1,mTrNr);
        for i = 1: size(avgSessionM,1)
            [mPeak(i), peak_index] = getPeakPostStim(avgSessionM(i,:), 15, 1, 2, 1);
            mBaseline(i) = mean(avgSessionM(i,1:30),2);
        end

        figure
        p1 = plot (1:mTrNr, mPeak, 'ko-', 'lineWidth', 1.5); hold on
        p2 = plot(1:mTrNr, mBaseline, 'k:', 'lineWidth', 1.5);
        [r,p,xfit,yfit] = fit_line(1:5, mPeak(1:5)');
        p3 =plot(xfit, yfit, 'r', 'LineWidth', 2);
        legend ([p1 p2], {'postStim', 'baseline'});
        legend boxoff

        xlabel('trial number'); ylabel('dF/F_0')
        title('monoc'); axis square; box off
        xlim([0 mTrNr]); ylim([-0.02 0.08]);
        set(findobj(gcf, 'type', 'axes'), 'FontSize', 12)
        savefig(gcf, 'peak_vs_MonocHitTrialNumber')
        saveas(gcf, 'peak_vs_MonocHitTrialNumber.png')
    end
    clear avgTempTrace indivTempTrace avgSessionB avgSessionM bTrNr mTrNr
end 

function [r,p,xfit,yfit] = fit_line(xdata,ydata)

coefficients = polyfit(xdata, ydata,1);
npoints = length(xdata); 
xfit = linspace(min(xdata), max(xdata), npoints);
yfit = polyval(coefficients, xfit);
%plot(xfit, yfit, '--k'); hold on;
[r,p] = corrcoef(xdata, ydata);
end