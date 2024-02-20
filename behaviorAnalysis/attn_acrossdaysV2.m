clear; close all; clc

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK\project'))

animal = 'M230717_1';
% days = {'09-Oct-2023' '10-Oct-2023' '12-Oct-2023' '13-Oct-2023' '15-Oct-2023_1' '17-Oct-2023' '19-Oct-2023'...
%     '20-Oct-2023' '25-Oct-2023' '26-Oct-2023'}; % this is WFI experiment ID format
sPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' animal '\behavior\attention_over_days'];
cd (sPath)
load ("key_attn.mat")

% daysB = {'2023-10-09', '2023-10-13', '2023-10-19', '2023-10-20'};
idx = cellfun(@(x) ~isempty(x), key.expInfoNr);
daysB = key.expID(idx(:,1));
daysM = key.expID(idx(:,2));

% attnIdxB = {'binoc', 'monoc', 'monoc', 'binoc', 'monoc', 'monoc', 'binoc', 'binoc', 'monoc', 'monoc'};
sessB =  key.expInfoNr(idx(:,1),1)';
sessM =  key.expInfoNr(idx(:,2),2)';

attn_acrossdays(animal, daysB, sessB, 'binoc', sPath);
attn_acrossdays(animal, daysM, sessM, 'monoc', sPath);

%%
function attn_acrossdays(animal, days, sess, flag, savePath)

keynew = struct;

for d = 1:length(days)
    keynew.(animal)(d).Date = datestr(days{d}(1:11), 'yyyy-mm-dd');
    [out, beh_all, ~, trTable, attn] = attention_effects_hitAndTrialAligned_EK_V2(animal,datestr(days{d}(1:11), 'yyyy-mm-dd'),  split(num2str(sess{d}))');
    close all
    keynew.(animal)(d).Out = beh_all;
    keynew.(animal)(d).attn = attn;
end
key = keynew;
%% Reaction Time aligned to hit trials
bRT = [];
mRT = [];

for d = 1:length(key.(animal))
    % for s = 1:length(key.(animal)(d).Sessions)
    % if ~isempty(key.(animal)(d).Sessions(s).Out)
    switch flag
        case  'binoc'
            bRT = [bRT; key.(animal)(d).attn.alignToBlockMeanRtB'];
            avg_RT_days = nanmean(bRT,1);
            std_RT_days = nanstd(bRT,0,1);

        case  'monoc'
            mRT = [mRT; key.(animal)(d).attn.alignToBlockMeanRtM'];
            avg_RT_days = nanmean(mRT,1);
            std_RT_days = nanstd(mRT,0,1);
    end
end 
    
%% Hit Rate aligned to hit trials
mdp =[];
bdp = [];

for d = 1:length(key.(animal))
    switch flag
        case 'binoc'
            bdp = [bdp; key.(animal)(d).attn.alignToBlockMeanHrB'];
            avg_dp_days = nanmean(bdp,1);
            std_dp_days = nanstd(bdp,0,1);

        case 'monoc'
            mdp = [mdp; key.(animal)(d).attn.alignToBlockMeanHrM];
            avg_dp_days = nanmean(mdp,1);
            std_dp_days = nanstd(mdp,0,1);

    end
end

%%
figure;
subplot (1,2,1)
errorbar(avg_RT_days, std_RT_days, '-o', 'LineWidth', 1.5, 'Color', [128 128 128]./255, 'CapSize', 0); hold on;
[r,p,xfit,yfit] = fit_line(1:5, avg_RT_days(1:5)');
plot(xfit, yfit, 'r', 'LineWidth', 2); 
ylabel(['Reaction Time (s)']);
xlabel(['Trial Number']);
ylim([0.25 0.7])
set(gca,'TickDir','out');
box off;
xlim([0 10]);
xticks([2:2:10]);
axis square

subplot (1,2,2)
errorbar(avg_dp_days, std_dp_days, '-o', 'LineWidth', 1.5, 'Color', [128 128 128]./255, 'CapSize', 0); hold on;
[r,p,xfit,yfit] = fit_line(1:5, avg_dp_days(1:5)');
plot(xfit, yfit, 'r', 'LineWidth', 2); 
ylabel(['Hit rate']);
xlabel(['Trial Number']);
set(gca,'TickDir','out');
box off;
axis square
ylim([0 1.2]);
xlim([0 10]);
xticks([2:2:10]);

sgtitle(flag)

cd(savePath)
savefig(gcf, ['behData_attn_' flag])


% %%
% subplot(3,2,4); hold on;
% % if b_slope_d > 0
% %     plot(avg_bdp_days, '-o', 'LineWidth', 1.2, 'Color', [128 128 128]./255); hold on;
% % elseif b_slope_d < 0
%     plot(avg_bdp_days, 'LineWidth', 2, 'Color', 'k'); hold on;
%       plot ([1:10], bdp(:,1:10), 'Color', [0.5 0.5 0.5])
% 
% % end
% [r,p,xfit,yfit] = fit_line(1:5, avg_bdp_days(1:5)');
% plot(xfit, yfit, 'b', 'LineWidth', 2); hold on;
% ylabel(['d prime']);
% xlabel(['Trial Number']);
% set(gca,'TickDir','out');
% box off;
% axis square
% ylim([0 3]);
% xlim([0 11]);
% xticks([2:2:10]);
% 
% 
% mhit =[];
% bhit= [];
% mmiss =[];
% bmiss= [];
% mfa =[];
% bfa= [];
% 
% for d = 1:length(key.(animal))
%     % for s = 1:length(key.(animal)(d).Sessions)
%     %     if ~isempty(key.(animal)(d).Sessions(s).Out)
%             mhit = [mhit; key.(animal)(d).Attn.m_firstHits'];
%             bhit = [bhit; key.(animal)(d).Attn.b_firstHits'];
%              mmiss = [mmiss; key.(animal)(d).Attn.m_firstMiss'];
%             bmiss = [bmiss; key.(animal)(d).Attn.b_firstMiss'];
%              mfa = [mfa; key.(animal)(d).Attn.m_firstfa'];
%             bfa = [bfa; key.(animal)(d).Attn.b_firstfa'];
%         % end
%     % end
%     % avg_mdp(d,:) = nanmean(mdp,1);
%     % std_mdp(d,:) = nanstd(mdp,0,1);
%     % avg_bdp(d,:) = nanmean(bdp,1);
%     % std_bdp(d,:) = nanstd(bdp,0,1);
%     % mdp = [];
%     % bdp = [];
% 
% end
% 
% m_firstHits = nanmean(mhit,1);
% m_firstMiss =nanmean(mmiss,1);
% m_firstfa = nanmean(mfa,1);
% b_firstHits = nanmean(bhit,1);
% b_firstMiss =nanmean(bmiss,1);
% b_firstfa = nanmean(bfa,1);
% 
% subplot(3,2,5); hold on;
% plot(m_firstHits, '-o', 'LineWidth', 1.2, 'Color', 'g'); hold on;
% plot(m_firstMiss, '-o', 'LineWidth', 1.2, 'Color', 'r'); hold on;
% plot(m_firstfa, '-o', 'LineWidth', 1.2, 'Color', 'c'); hold on;
% ylabel(['rate (%)']);
% xlim([0 11]);
% xticks([2:2:10]);
% set(gca,'TickDir','out');
% xlabel(['Trial Number']);
% axis square
% 
% subplot(3,2,6); hold on;
% plot(b_firstHits, '-o', 'LineWidth', 1.2, 'Color', 'g'); hold on;
% plot(b_firstMiss, '-o', 'LineWidth', 1.2, 'Color', 'r'); hold on;
% plot(b_firstfa, '-o', 'LineWidth', 1.2, 'Color', 'c'); hold on;
% ylabel(['rate (prop.)']);
% xlim([0 11]);
% xticks([2:2:10]);
% set(gca,'TickDir','out');
% xlabel(['Trial Number']);
% 
% axis square
end 