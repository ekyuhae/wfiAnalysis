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

Animal = 'M230717_1';
hva = 'LM_PM';
sPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior\attention_over_days'];
cd (sPath)
load ("key_attn.mat")
if ~isempty(hva)
    sPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior\attention_over_days_' hva];
    isfolder(sPath);
    mkdir(sPath);
end

expID = key.expID;
flag = 'matched'; % whether retinotopically matched ('matched'), ROI switched ('roiUnmatched'), or block switched ('blockUnmatched')
% attIdx as a session that had attentional effects (binoc or monoc block)
% attIdx = {'b','m','m','b','m','m','b','b','m','m'}; % M230831
attIdx = {'b','b','b','m','b','b','m','b','b','b','m'}; % M230717

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
data = struct;
binocTemp = []; monocTemp =[];
for d = 1:length(expID)
    if isempty(hva)
        cd([fPath filesep Animal filesep 'behavior' filesep expID{d} '\attention_over_sessions'])
    else 
         cd([fPath filesep Animal filesep 'behavior' filesep expID{d} '\attention_over_sessions_' hva])
    end 
    switch flag
        case 'matched'
            load('hit_aligned_dff.mat', 'avgTempTrace');
        case 'roiUnmatched'
            load('hit_aligned_dff_ROIswitched.mat', 'avgTempTrace');
        case 'blockUnmatched'
            load('hit_aligned_dff_blockSwitched.mat', 'avgTempTrace');
    end
    if ~strcmp(flag, 'blockUnmatched')
        if isfield(avgTempTrace, "binoc") && attIdx{d} == 'b'
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

        if isfield(avgTempTrace, "monoc") && attIdx{d} == 'm'
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
    else  % for block unmatched case; e.g. bino attn session and looking at response in monoc block in that session
        if isfield(avgTempTrace, "binoc") && attIdx{d} == 'm' 
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

        if isfield(avgTempTrace, "monoc") && attIdx{d} == 'b'
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
fr = 30/2;
t = 0:1/fr:90/fr; % convert frame into timescale
t = t - t(fr*2+1);
cd (sPath)
    figure ;
    if height(avgDaysB) < 10
        k = height(avgDaysB);
    else
        k = 10;
    end
    for i = 1: k
        subplot(1,k,i)
        plot(t, avgDaysB(i,:), 'LineWidth', 1.5, 'color','k');
        xline(t(fr*2+1), ':k');
        title(['hit no.' num2str(i)]);
        xlabel('time (s)'); ylabel('dF/F_0')
        axis square; box off
        ylim([-0.05 0.08])
    end
    sgtitle(['Binoc Attentional block across days']) % MONOC ATTENTION SESSION FOR BLOCK UNMATCHED

    savefig(gcf, ['dffAcrossBinocHitTrials_' flag] )
    saveas(gcf, ['dffAcrossBinocHitTrials_' flag '.png'])



    figure;
    if height(avgDaysM) < 10
        k = height(avgDaysM);
    else
        k = 10;
    end
    for i = 1: k
        subplot(1,k,i)
        plot(t, avgDaysM(i,:), 'LineWidth', 1.5, 'color','k');

        xline(t(fr*2+1), ':k');
        title(['monoc hit no.' num2str(i)]);
        xlabel('time (s)'); ylabel('dF/F_0')
        axis square; box off
        ylim([-0.05 0.08])
    end
    sgtitle(['Monoc Attentional block across days'])
    savefig(gcf, ['dffAcrossMonocHitTrials_' flag])
    saveas(gcf, ['dffAcrossMonocHitTrials_' flag '.png'])

%%
bTrNr = 15; mTrNr = 15;
if exist("bTrNr", "var")
    bTrNr = max(bTrNr); bPeak = nan(1,bTrNr); bPeak_std = nan(1,bTrNr); bBaseline = nan(1,bTrNr); bBaselineStd = nan(1,bTrNr);
    for i = 1: size(avgDaysB,1)
        [bPeak(i), peak_index] = getPeakPostStim(avgDaysB(i,:), 15, 1, 2, 1);
        [bPeak_std(i), peak_index] = getPeakPostStim(avgDaysB_std(i,:), 15, 1, 2, 1);
        % [bBaseline(i), ~] = getPeakPostStim(avgDaysB(i,:), 15, 2, 0, 1);
        % [bBaselineStd(i), ~] = getPeakPostStim(avgDaysB_std(i,:), 15, 2, 0, 1);

        bBaseline(i) = mean(avgDaysB(i,1:30),2);
        bBaselineStd(i) = std(avgDaysB(i,1:30),0,2);
    end

    figure;

    % p1 = plot(1:bTrNr, bPeak, 'ko-'); hold on
    p1 = errorbar (1:bTrNr, bPeak(1:bTrNr), bPeak_std(1:bTrNr), 'ko-', 'CapSize', 0); hold on

    p2 = errorbar(1:bTrNr, bBaseline(1:bTrNr), bBaselineStd(1:bTrNr), 'k:', 'CapSize', 0);
    [r,p,xfit,yfit] = fit_line(1:5, bPeak(1:5)');
    p3 = plot(xfit, yfit, 'r', 'LineWidth', 2);
    legend ([p1 p2], {'postStim', 'baseline'})

    % errorbar (1:bTrNr, bBaseline, bBaselineStd ,'-o', 'CapSize',0, 'color', [0.5 0.5 0.5])
    xlabel('trial number'); ylabel('dF/F_0')
    title('binoc'); axis square; box off
    xlim([0 bTrNr]); ylim([-0.03 0.07])
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
        % [mBaseline(i), ~] = getPeakPostStim(avgDaysM(i,:), 15, 2, 0, 1);
        % [mBaselineStd(i), ~] = getPeakPostStim(avgDaysM_std(i,:), 15, 2, 0, 1);

        mBaseline(i) = mean(avgDaysM(i,1:30),2);
        mBaselineStd(i) = std(avgDaysM(i,1:30),0,2);

    end

    figure
    p1 = errorbar (1:mTrNr, mPeak(1:mTrNr), mPeak_std(1:mTrNr), 'ko-', 'CapSize', 0); hold on
    p2 = errorbar(1:mTrNr, mBaseline(1:mTrNr), mBaselineStd(1:mTrNr), 'k:', 'CapSize', 0);
    [r,p,xfit,yfit] = fit_line(1:5, mPeak(1:5)');
    p3 =plot(xfit, yfit, 'r', 'LineWidth', 2);
    legend ([p1 p2], {'postStim', 'baseline'})

    xlabel('trial number'); ylabel('dF/F_0')
    title('monoc'); axis square; box off
    xlim([0 mTrNr]); ylim([-0.03 0.07])
    if strcmp(flag, 'matched')
        savefig(gcf, 'peak_vs_MonocHitTrialNumber')
        saveas(gcf, 'peak_vs_MonocHitTrialNumber.png')
    else
        savefig(gcf, ['peak_vs_MonocHitTrialNumber_' flag])
        saveas(gcf, ['peak_vs_MonocHitTrialNumber_' flag '.png'])
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