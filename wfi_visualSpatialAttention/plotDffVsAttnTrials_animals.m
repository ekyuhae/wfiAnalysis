%% WFI Behavior Experiment Pipeline - Analysis for ATTENTIONAL EFFECTS
% analysis after running wfiBehavior (analysis for visual detection; df/f
% PSTHs for different tial types) since it requires
% 'analyed_active_(sessNr).mat' for pulling ROI
% pull the data from the session that has attentional signs in binoc and/or
% monoc blocks
%
% EK Dec 23

clc; clear; close all

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

Animal = {'M230717_1' 'M230831_1'};
sPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior\attention_over_days'];
cd (sPath)
load ("key_attn.mat")

% expID = {'09-Oct-2023' '10-Oct-2023' '12-Oct-2023' '13-Oct-2023' '15-Oct-2023_1' '17-Oct-2023' '19-Oct-2023'...
%     '20-Oct-2023' '23-Oct-2023' '25-Oct-2023' '26-Oct-2023'}; % this is WFI experiment ID format
expID = key.expID;
flag = 'blockUnmatched'; % whether retinotopically matched ('matched'), ROI switched ('roiUnmatched'), or block switched ('blockUnmatched')
% sPath = 'Y:\haider\Data\analyzedData\EKK\WFI\M230831_1\behavior\attention_over_days';

% expID = {'15-Oct-2023_1' '20-Oct-2023' '23-Oct-2023' '25-Oct-2023' '26-Oct-2023'}; % this is WFI experiment ID format

% hva = []; %higher visual areas being analyzed (i.e. 'RL_AM' or 'LM_PM' or '' for V1)
% % a.binoc ={[2]; [5 6]; [1 5 6]; [5]; [2 4 5]; [2 8]; [6 7 9]; []; [5 7]; [9]; [5]; []; []};
% a.binoc ={[3 4 5]; [9]; [5]; []; []};
%
% % a.monoc = {[]; []; [1 2]; [1 4]; [1 2 5 6]; [1 3]; [7]; [1 4 6 7]; [1 3]; [4 7 9]; [3 6]; [1 3 4 6]; [1 2 5 6]} ; % enter the imaging session number(X) from Frames_2_xxx_xxx_uint16_000X
% a.monoc = {[4]; [4 7 9]; [3 6]; [1 3 4 6]; [1 2 5 6]} ; % enter the imaging session number(X) from Frames_2_xxx_xxx_uint16_000X
%
% b.binoc = {[6 7 9]; [9]; [5]; []; []};
% b.monoc = {[7]; [4 7 9]; [3 6]; [1 3 4 6]; [1 2 5 6]} ;
% sessNr = [a.binoc, a.monoc];
% expInfoNr = [b.binoc, b.monoc];
% clear a b
%% pulling the analyzed mat file (output from running wfiBehavior) and get df/f averages across sessions for each day
% avgSessionB/M; average across sessions for each day (1xdays)
fPath ='Y:\haider\Data\analyzedData\EKK\WFI\';
for iAnimal = 1:length(Animal)
data = struct;
binocTemp = []; monocTemp =[];
for d = 1:length(expID)
    cd([fPath filesep Animal filesep 'behavior' filesep expID{d} '\attention_over_sessions'])
    switch flag
        case 'matched'
            load('hit_aligned_dff.mat', 'avgTempTrace');
        case 'roiUnmatched'
            load('hit_aligned_dff_ROIswitched.mat', 'avgTempTrace');
        case 'blockUnmatched'
            load('hit_aligned_dff_blockSwitched.mat', 'avgTempTrace');
    end
    if isfield(avgTempTrace, "binoc")
        nrAlignedHits = max(arrayfun(@(x) size(x.binoc,1), avgTempTrace));
        for trial= 1:nrAlignedHits
            for i = 1:length(avgTempTrace)

                if trial <= height(avgTempTrace(i).binoc)
                    binocTemp = [binocTemp; avgTempTrace(i).binoc(trial,:)];
                else
                    binocTemp = binocTemp;
                end
            end
            avgSessionB{d}(trial,:) = mean(binocTemp,1);
            avgSessionB_std{d}(trial,:) = std(binocTemp,0,1);

            binocTemp =[];
        end
    end

    if isfield(avgTempTrace, "monoc")
        nrAlignedHits = max(arrayfun(@(x) size(x.monoc,1), avgTempTrace));
        for trial= 1:nrAlignedHits
            for i = 1:length(avgTempTrace)

                if trial <= height(avgTempTrace(i).monoc)
                    monocTemp = [monocTemp; avgTempTrace(i).monoc(trial,:)];
                else
                    monocTemp = monocTemp;
                end
            end
            avgSessionM{d}(trial,:) = mean(monocTemp,1);
            avgSessionM_std{d}(trial,:) = std(monocTemp,0,1);

            monocTemp =[];
        end
    end
end

%%

binocTemp = []; monocTemp =[];

nrAlignedHits = max(cellfun(@(x) height(x), avgSessionB));
for trial= 1:nrAlignedHits
    for i = 1:length(avgSessionB)
        if trial <= height(avgSessionB{i})
            binocTemp = [binocTemp; avgSessionB{i}(trial,:)];
        else
            binocTemp = binocTemp;
        end
    end
    avgDaysB(trial,:) = mean(binocTemp,1);
    avgDaysB_std(trial,:) = std(binocTemp,0,1);

    binocTemp =[];
end

nrAlignedHits = max(cellfun(@(x) height(x), avgSessionM));
for trial= 1:nrAlignedHits
    for i = 1:length(avgSessionM)
        if  trial <= height(avgSessionM{i})
            monocTemp = [monocTemp; avgSessionM{i}(trial,:)];
        else
            monocTemp = monocTemp;
        end
    end
    avgDaysM(trial,:) = mean(monocTemp,1);
    avgDaysM_std(trial,:) = std(monocTemp,0,1);

    monocTemp =[];
end

%%
% fr = 30/2;
% t = 0:1/fr:90/fr; % convert frame into timescale
% t = t - t(fr*2+1);
% cd (sPath)
%     figure ;
%     if height(avgDaysB) < 10
%         k = height(avgDaysB);
%     else
%         k = 10;
%     end
%     for i = 1: k
%         subplot(1,k,i)
%         plot(t, avgDaysB(i,:), 'LineWidth', 1.5, 'color','k');
%         xline(t(fr*2+1), ':k');
%         title(['hit no.' num2str(i)]);
%         xlabel('time (s)'); ylabel('dF/F_0')
%         axis square; box off
%         ylim([-0.05 0.08])
%     end
%     % sgtitle(['Binoc Attentional block: sess ' strsplit(num2str(sessNr{cDay}), ' ')])
% 
%     % savefig(gcf, 'dffAcrossBinocHitTrials')
%     % saveas(gcf, 'dffAcrossBinocHitTrials.png')
% 
% 
% 
%     figure;
%     if height(avgDaysM) < 10
%         k = height(avgDaysM);
%     else
%         k = 10;
%     end
%     for i = 1: k
%         subplot(1,k,i)
%         plot(t, avgDaysM(i,:), 'LineWidth', 1.5, 'color','k');
% 
%         xline(t(fr*2+1), ':k');
%         title(['monoc hit no.' num2str(i)]);
%         xlabel('time (s)'); ylabel('dF/F_0')
%         axis square; box off
%         ylim([-0.05 0.08])
%     end
%     % sgtitle(['Monoc Attentional block: sess ' strsplit(num2str(sessNr{cDay}), ' ')])
%     % savefig(gcf, 'dffAcrossMonocHitTrials')
%     % saveas(gcf, 'dffAcrossMonocHitTrials.png')

%%
bTrNr = 15; mTrNr = 15;
if exist("bTrNr", "var")
    bTrNr = max(bTrNr); bPeak = nan(1,bTrNr); bPeak_std = nan(1,bTrNr); bBaseline = nan(1,bTrNr); bBaselineStd = nan(1,bTrNr);
    for i = 1: size(avgDaysB,1)
        [bPeak(i), peak_index] = getPeakPostStim(avgDaysB(i,:), 15, 1, 2, 1);
        [bPeak_std(i), peak_index] = getPeakPostStim(avgDaysB_std(i,:), 15, 1, 2, 1);
        bBaseline(i) = mean(avgDaysB(i,1:30),2);
        % bBaselineStd(i) = std(avgSessionB(i,1:30),0,2);
    end
%% across animals 
tempB(:,:,iAnimal) = avgDaysB;
tempM(:,:,iAnimal) = avgDaysM;
avgAniB = mean(tempB,3); 
avgAniM = mean(tempM,3);

    figure;

    % p1 = plot(1:bTrNr, bPeak, 'ko-'); hold on
    p1 = errorbar (1:bTrNr, bPeak(1:bTrNr), bPeak_std(1:bTrNr), 'ko-', 'CapSize', 0); hold on

    p2 = plot(1:bTrNr, bBaseline(1:bTrNr), 'k:');
    [r,p,xfit,yfit] = fit_line(1:5, bPeak(1:5)');
    p3 = plot(xfit, yfit, 'r', 'LineWidth', 2);
    legend ([p1 p2], {'postStim', 'baseline'})

    % errorbar (1:bTrNr, bBaseline, bBaselineStd ,'-o', 'CapSize',0, 'color', [0.5 0.5 0.5])
    xlabel('trial number'); ylabel('dF/F_0')
    title('binoc'); axis square; box off
    xlim([0 bTrNr]); ylim([-0.04 0.08])
    if strcmp(flag, 'matched')
        savefig(gcf, 'peak_vs_BinocHitTrialNumber')
        saveas(gcf, 'peak_vs_BinocHitTrialNumber.png')
    else
        savefig(gcf, ['peak_vs_BinocHitTrialNumber_' flag])
        saveas(gcf, ['peak_vs_BinocHitTrialNumber_' flag '.png'])
    end

end

if exist("mTrNr", "var")
    mTrNr = max(mTrNr); mPeak = nan(1,mTrNr);mPeak_std = nan(1,mTrNr); mBaseline = nan(1,mTrNr); mBaselineStd = nan(1,mTrNr);
    for i = 1: size(avgDaysM,1)
        [mPeak(i), peak_index] = getPeakPostStim(avgDaysM(i,:), 15, 1, 2, 1);
        [mPeak_std(i), peak_index] = getPeakPostStim(avgDaysM_std(i,:), 15, 1, 2, 1);

        mBaseline(i) = mean(avgDaysM(i,1:30),2);
    end

    figure
    p1 = errorbar (1:mTrNr, mPeak(1:mTrNr), mPeak_std(1:mTrNr), 'ko-', 'CapSize', 0); hold on
    p2 = plot(1:mTrNr, mBaseline(1:mTrNr), 'k:');
    [r,p,xfit,yfit] = fit_line(1:5, mPeak(1:5)');
    p3 =plot(xfit, yfit, 'r', 'LineWidth', 2);
    legend ([p1 p2], {'postStim', 'baseline'})

    xlabel('trial number'); ylabel('dF/F_0')
    title('monoc'); axis square; box off
    xlim([0 mTrNr]); ylim([-0.04 0.08])
    if strcmp(flag, 'matched')
        savefig(gcf, 'peak_vs_MonocHitTrialNumber')
        saveas(gcf, 'peak_vs_MonocHitTrialNumber.png')
    else
        savefig(gcf, ['peak_vs_MonocHitTrialNumber_' flag])
        saveas(gcf, ['peak_vs_MonocHitTrialNumber_' flag '.png'])
    end

end
end 
function [r,p,xfit,yfit] = fit_line(xdata,ydata)

coefficients = polyfit(xdata, ydata,1);
npoints = length(xdata);
xfit = linspace(min(xdata), max(xdata), npoints);
yfit = polyval(coefficients, xfit);
%plot(xfit, yfit, '--k'); hold on;
[r,p] = corrcoef(xdata, ydata);
end