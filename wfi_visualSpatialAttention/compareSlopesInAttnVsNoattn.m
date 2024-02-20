clear; close all; clc

addpath(genpath("Y:\haider\Code\behaviorAnalysis"));
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
Animal = 'M230717_1';
hva = {'', 'LM_PM', 'RL_AM'};
sPath = ['Y:\haider\Data\analyzedData\EKK\WFI\' Animal '\behavior\attention_over_days'];
cd (sPath)
load ("key_attn.mat")
% attnblk = 'B';

expID = key.expID;
flag = 'matched'; % whether retinotopically matched ('matched'), ROI switched ('roiUnmatched'), or block switched ('blockUnmatched')
% attIdx as a session that had attentional effects (binoc or monoc block)
% attIdx = {'b','m','m','b','m','m','b','b','m','m'}; % M230831
attIdx = {'b','b','b','m','b','b','m','b','b','b','m'}; % M230717


%%
for area = 1:length(hva)
    for d = 1:length(expID)
        fPath ='Y:\haider\Data\analyzedData\EKK\WFI\';
        % data = struct;
        % binocTemp = []; monocTemp =[];

        if isempty(hva{area})
            cd([fPath filesep Animal filesep 'behavior' filesep expID{d} '\attention_over_sessions'])
        else
            cd([fPath filesep Animal filesep 'behavior' filesep expID{d} '\attention_over_sessions_' hva{area}])
        end
        if attIdx{d} == 'b'
            fig1 = openfig(['peak_vs_BinocHitTrialNumber.fig']);
            if exist('peak_vs_MonocHitTrialNumber_blockSwitched.fig')
                fig2 = openfig('peak_vs_MonocHitTrialNumber_blockSwitched.fig');
            end
        else
            fig1 = openfig(['peak_vs_MonocHitTrialNumber.fig']);
            if exist('peak_vs_BinocHitTrialNumber_blockSwitched.fig')
                fig2 = openfig('peak_vs_BinocHitTrialNumber_blockSwitched.fig');
            end

        end
        if ~isvarname(fig1) && ~isvarname(fig2)
            axes1 = findobj(fig1, 'Type', 'axes');
            data1 = get(axes1, 'Children'); % Get plot objects
            x_data1 = get(data1, 'XData'); % Get x data
            y_data1 = get(data1, 'YData'); % Get y data

            % Figure 2
            axes2 = findobj(fig2, 'Type', 'axes');
            data2 = get(axes2, 'Children'); % Get plot objects
            x_data2 = get(data2, 'XData'); % Get x data
            y_data2 = get(data2, 'YData'); % Get y data

            % calculate slopes
            slp_matched{area}(d) = y_data1{1}(2)-y_data1{1}(1);
            slp_unmatched{area}(d) =y_data2{1}(2)-y_data2{1}(1);
        end

        clear data1 data2 x_data1 x_data2 y_data1 y_data2
    end
end

%%
close all
% Plot the extracted data points from both figures together
figure;
for area = 1:length(hva)
    for d = 1:length(slp_unmatched{area})
        % if area == 1 % V1
            plot([1 2], [slp_matched{area}(d) slp_unmatched{area}(d)], 'o-','color',[0.76,0.87,0.78] ,'lineWidth', 1,'MarkerSize',5, 'MarkerFaceColor',[0.16,0.38,0.27])
        % elseif area ==2 % HVA1
        %     plot([1 2], [slp_matched{area}(d) slp_unmatched{area}(d)], 'o-', 'color', [0.7 0.1 0.3], 'MarkerFaceColor',[0.5 0.5 0.5])
        % else % HVA2
        %     plot([1 2], [slp_matched{area}(d) slp_unmatched{area}(d)], 'o-', 'color', [0.3 0.7 0.2], 'MarkerFaceColor',[0.5 0.5 0.5])
        % end
        hold on;
    end
end
box off; axis square; setfig(); xlim([0.5 2.5]); ylim([-0.008 0.01])
ylabel('Slope')
x_labels = {'Attention Block' 'Non-Attention Block'};
set(gca, 'XTick', 1:length(x_labels), 'XTickLabel', x_labels);
%%
% p_values = zeros(1, length(m_p));
% for i = 2:length(m_p)
%     [~, p_values(i)] = ttest(single(aa{2}(i,:))', a{2}(i,:)');
% end
% 
% % Add p-values as text annotations
% p_value_threshold = 0.05;
% for i = 2:5
%     if ~isnan(p_values(i)) && p_values(i) < p_value_threshold
%         text(x(2) + 0.1, m_a(i,1)*100, "*", 'FontSize', 12); hold on;
%     end
% end
% text(x(1), m_p + 0.005, num2str(p_values', '%0.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'r');
% text(x(2), m_a + 0.005, num2str(p_values', '%0.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'b');

cd([fPath filesep Animal filesep 'behavior']);
savefig('slopeComparison_attnBlockVsNonattnBlock.fig');

 