close all; clear all; clc;
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\';
Animal = {'M221118_1','M230130_2','M230131_1'};
cont = [100 75 50 25 5 0];
full = [1 0 0 0 0 1]; five = [0 0 0 0 1 1]; multi = [1 1 1 1 0 1];

data = struct();
data.(Animal{1})(1,:) = {'02-Feb-2023_1', '03-Feb-2023', '08-Feb-2023_2', '15-Feb-2023_2', '18-Feb-2023', '09-Mar-2023'};
data.(Animal{1})(2,:) ={full, full, full, multi, multi, five};
data.(Animal{2})(1,:) = {'20-Feb-2023', '24-Feb-2023', '09-Mar-2023'};
data.(Animal{2})(2,:) ={full, multi, five};
data.(Animal{3})(1,:) = {'01-Mar-2023', '01-Mar-2023_1'};
data.(Animal{3})(2,:) ={full, multi};

combinedROIs = true; % picking [best-1 best best+1] ROI response

%%
traces = struct();
for i = 1:length(Animal)
   for exp= 1:numel(data.(Animal{i})(1,:))
    temp  = compileIndividTrial(Animal{i}, data.(Animal{i}){1,exp}, fPath, sum(data.(Animal{i}){2,exp}));
    % dFFtraces = [dFFtraces; temp];
    traces(i,exp).trace = temp;
   end
end
%%
dat = cell(length(cont), 17);
dFFtraces =cell(length(cont), 17);
for i = 1:length(Animal)
    % locations
    for exp= 1:numel(data.(Animal{i})(1,:))
        for jj = 1: length(traces(1,1).trace)
            traces(i,exp).trace{end,jj} = traces(i,exp).trace{end,1}; % for 0 contrast
        end

        if isequal(data.(Animal{i}){2,exp},full)

            dat(1,:) = traces(i,exp).trace(1,:);
            dat(end,:) = traces(i,exp).trace(end,:);
        elseif isequal(data.(Animal{i}){2,exp},multi)
            dat(1:4,:) = traces(i,exp).trace(1:4,:);
            dat(6,:) = traces(i,exp).trace(end,:);
        elseif isequal(data.(Animal{i}){2,exp},five)
            dat(5,:) = traces(i,exp).trace(1,:);
            dat(6,:) = traces(i,exp).trace(end,:);
        end

        dFFtraces = cellfun(@(x,y) [x;y], dFFtraces, dat, 'UniformOutput',false);
        dat = cell(length(cont), 17);
    end
end
%%
if combinedROIs
       trace5 = [];
    trace25 = [];
    trace50 = [];
    trace75 = [];
    trace100 = [];
    % [s1, s2] = size(dFFtraces);
    % for j= 1: s1
    %     for jj = 1:s2
    %         if ~isempty(dFFtraces(j,jj).dFFtrace)
                for i = 1:17
                    if i ==1
                        trace5 = [trace5; dFFtraces{5,i}; dFFtraces{5,i+1} ];
                        trace25 = [trace25; dFFtraces{4,i}; dFFtraces{4,i+1}];
                        trace50 = [trace50; dFFtraces{3,i}; dFFtraces{3,i+1}];
                        trace75 = [trace75; dFFtraces{2,i}; dFFtraces{2,i+1}];
                        trace100 = [trace100; dFFtraces{1,i}; dFFtraces{1,i+1}];
                    elseif i ==17
                        trace5 = [trace5; dFFtraces{5,i-1}; dFFtraces{5,i}];
                         trace25 = [trace25; dFFtraces{4,i-1}; dFFtraces{4,i}];
                          trace50 = [trace50; dFFtraces{3,i-1}; dFFtraces{3,i}];
                           trace75 = [trace75; dFFtraces{2,i-1}; dFFtraces{2,i}];
                            trace100 = [trace100; dFFtraces{1,i-1}; dFFtraces{1,i}];
                    else
                        trace5 = [trace5; dFFtraces{5,i-1}; dFFtraces{5,i};dFFtraces{5,i+1}];
                        trace25 = [trace25; dFFtraces{4,i-1}; dFFtraces{4,i};dFFtraces{4,i+1}];
                        trace50 = [trace50; dFFtraces{3,i-1}; dFFtraces{3,i};dFFtraces{3,i+1}];
                        trace75 = [trace75; dFFtraces{2,i-1}; dFFtraces{2,i};dFFtraces{2,i+1}];
                        trace100 = [trace100; dFFtraces{1,i-1}; dFFtraces{1,i};dFFtraces{1,i+1}];
                    end
                end

     tracedata{5} = trace5; tracedata{4} = trace25; tracedata{3} = trace50; tracedata{2} = trace75; tracedata{1} = trace100;   
    noisetrace = dFFtraces{6,1};    

    % noisetrace =[];
    % for j = 1:s1
    %     for jj = 1:s2
    %         if ~isempty(dFFtraces(j,jj).dFFtrace)
    %             noisetrace = [noisetrace; dFFtraces(j,jj).dFFtrace{2,1}];
    % 
    %         end
    %     end
    % end

else % original ROI traces

    trace5 = [];
    trace25 = [];
    trace50 = [];
    trace75 = [];
    trace100 = [];
    % [s1, s2] = size(dFFtraces);
    % for j= 1: s1
    %     for jj = 1:s2
    %         if ~isempty(dFFtraces(j,jj).dFFtrace)
    for i = 1:17
        trace5 = [trace5; dFFtraces{5,i}];
        trace25 = [trace25; dFFtraces{4,i}];
        trace50 = [trace50; dFFtraces{3,i}];
        trace75 = [trace75; dFFtraces{2,i}];
        trace100 = [trace100; dFFtraces{1,i}];
    end
    tracedata{5} = trace5; tracedata{4} = trace25; tracedata{3} = trace50; tracedata{2} = trace75; tracedata{1} = trace100; 
    %         end
    %     end
    % end
noisetrace = dFFtraces{6,1};    
%     noisetrace =[];
%     for j = 1:s1
%         for jj = 1:s2
%             if ~isempty(dFFtraces(j,jj).dFFtrace)
%                 % noisetrace = [noisetrace; dFFtraces(j,jj).dFFtrace{2,1}];
% 
% 
%             end
%         end
%     end
end
 % a = mean(noisetrace,2);
 % figure; plot(1:length(a),a)

%% 5 perc contrast 
% % whole trial 
%    for i = 1:size (noisetrace,1)
%         noise(i) =getPeakPostStim(noisetrace(i,:), 15,1, 2);
%     end
%%
for k = 1:length(tracedata)
    for i = 1: size(tracedata{k},1)
        tracedata1{k}(i,:) = tracedata{k}(i,:) - min(tracedata{k}(i,:));
        [signal{k}(i), peak_index_5perc(i)] = getPeakPostStim(tracedata1{k}(i,:), 15,1, 2);
        noise_fromstd{k}(i) = std(tracedata1{k}(i,24:30),0,2);
    end
 
    % figure; histogram(signal{k},50,'Normalization', 'probability'); hold on; histogram(noise_fromstd{k}',50,'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha',0.5)
    % legend([num2str(cont(k)) '% contrast'] , '0% contrast'); box off; legend('boxoff'); axis square
    % title ([num2str(cont(k)) '% contrast (' num2str(length(signal{k})) ' trials)']); xlabel('\DeltaF/F_0'); ylabel ('probability')
end 
leg= cell(5,1);
figure;
for k = 1:length(tracedata)
% noise_vec = reshape(cattrace(:,15:29), [], 1);
% noise_vec = reshape(noisetrace(:,31:end), [], 1);
  
    signal_vec = signal{k}';
    noise_vec = noise_fromstd{k}';
    pred = [signal_vec; noise_vec];
    resp = logical([repmat(1,length(signal_vec),1); repmat(0,length(noise_vec), 1)]); %index indicating signal (1) and noise (0)
    mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit'); %creation of distribution model
    scores = mdl.Fitted.Probability;
    outcome1 = cell(length(signal_vec),1);
    outcome1(:) = {'signal'};
    outcome2 = cell(length(noise_vec),1);
    outcome2(:) = {'noise'};
    outcome = [outcome1; outcome2];
    [X, Y, T, AUC(k)] = perfcurve(outcome, scores, 'signal');
    plot(X,Y); hold on; 
    xlabel('False positive rate');
    ylabel('True positive rate'); label = [num2str(cont(k)) '% contrast (' num2str(length(signal{k})) ' trials, AUC = ' num2str(AUC(k)) ')'];
    leg{k} = label;
    SNR{k} = signal{k}./noise_fromstd{k};
end 
plot([0,1], [0,1], 'k--') 
legend(leg)

cd([fPath filesep 'stats' filesep 'barmapping'])
save(['SNRandROC_' num2str(combinedROIs) '.mat'],"SNR","AUC","tracedata")
%%
for k = 1:5
snr(k) = mean(SNR{k}); sd(k)= std(SNR{k}, 0,2);
end 
figure; errorbar(cont(1:5), snr, sd, 'ko:', 'LineWidth', 1 );
% figure;  fillOut = fill([cont(1:5) fliplr(cont(1:5))],[snr+sd fliplr(snr-sd)],'k', 'FaceAlpha', 0.1,'linestyle','none'); hold on;
% lineOut = plot(cont(1:5),snr,'Color', 'k','linewidth',1.5);
ylabel('SNR'); xlabel('Contrast level(%)'); axis square; box off
ylim([-1 10])
%%
AUC_3rois = AUC; 
load("SNRandROC_0.mat")
figure; plot(cont(1:5), AUC, 'ko:','LineWidth', 2); hold on; plot(cont(1:5), AUC_3rois, 'o:','LineWidth', 2);
ylabel('AUROC'); xlabel('Contrast level(%)'); axis square; box off; legend('Original ROI', 'Combined ROIs')
% ylim([0.55 0.9]);

%%
trial_counts = [15,20,30,40,50,60,70,90];
 % trial_counts = [30,40,50,60,70,80,100, 200, 500, 1000];
n_subsamples = length(trial_counts);

for k = 1:5 % contrast: 100 ~ 0 
for run = 1:5000
    % Loop over the number of subsamples
    subsampled_data = cell(n_subsamples, 1);
    for i = 1:n_subsamples
        % Use datasample to randomly subsample trial_counts(i) trials with replacement
        subsampled_data{i} = datasample(tracedata1{k}, trial_counts(i),1, 'Replace', true);
        % subsampled_data_noise{i} = subsampled_data{i}(:,15:29);
        % subsampled_data_noise_mean{i} = mean(subsampled_data_noise{i},1);
        % subsampled_data_noise{i} = getPeakPostStim(cattrace(i,:), 15, 0.5, 1.5);
        for s = 1: trial_counts(i)
            [subtrace{i}(s), peak_index_5perc_sub(s)] = getPeakPostStim(subsampled_data{i}(s,:), 15, 1, 2);
            % [subnoise{i}(s), peak_index_noise_sub(s)] = getPeakPostStim(subsampled_data{i}(s,:), 15, 0.5, 1.5);
            subnoise{i}(s)=std(subsampled_data{i}(s,24:30),0,2);
        end
        % for s = 1: trial_counts(i)
        %     % subfiveperc{i}(s) = max(subsampled_data{i}(s,31:35));
        %     subsampled_data_noise{i}(s) = max(subsampled_data{i}(s,28:30));
        % end
%         signal_vec = subtrace{i}';
%         noise_vec = subnoise{i}';
%         pred = [signal_vec; noise_vec]; %put signal and noise in single array
%         resp = logical([repmat(1,length(signal_vec),1); repmat(0,length(noise_vec), 1)]); %index indicating signal (1) and noise (0)
%         mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit'); %creation of distribution model
%         scores = mdl.Fitted.Probability;
%         outcome1 = cell(length(signal_vec),1);
%         outcome1(:) = {'signal'};
%         outcome2 = cell(length(noise_vec),1);
%         outcome2(:) = {'noise'};
%         outcome = [outcome1; outcome2];
%     
%         [roc{k}(i).X, roc{k}(i).Y, roc{k}(i).T, roc{k}(i).AUC] = perfcurve(outcome,pred,'signal');
        SNR_sub{i,run} = subtrace{i}./subnoise{i};
% plot(roc(i).X,roc(i).Y);
% xlabel('False positive rate');
% ylabel('True positive rate'); hold on;
%         AUC_runs{k}(i, run) = roc{k}(i).AUC;
    end 
    clear subsampled_data 
end 
SNR_sub_contrasts_mean{k} = cell2mat(cellfun(@(x) mean(x),SNR_sub, 'UniformOutput',false));
SNR_sub_contrasts_std{k} = cell2mat(cellfun(@(x) std(x,0,2),SNR_sub, 'UniformOutput',false));
end 

figure;
k = linspace(0,0.6, 5);
for i = 1:5

 color = [k(i), k(i), k(i)];
    means = mean(SNR_sub_contrasts_mean{i},2);
sds = std(SNR_sub_contrasts_mean{i},0,2);
% errorbar(trial_counts, means, sds, ':o', 'Color', color, 'LineWidth',1.5); hold on
 fillOut = fill([trial_counts fliplr(trial_counts)],[(means+sds)' fliplr((means-sds)')],color, 'FaceAlpha', 0.1,'linestyle','none'); hold on;
 lineOut(i) = plot(trial_counts,means,'Color',color,'linewidth',1.5);
end 
ylabel('SNR'); xlabel('Number of Trials'); axis square; box off;
legend(lineOut, {'100%', '75%', '50%', '25%', '5%'})
%%
figure;
for i = 1: length(AUC_runs)
auc_mean{i} = mean(AUC_runs{i}, 2);
auc_std{i} = std(AUC_runs{i},1,2);
% auc_min = min(AUC_runs,[],2);
% auc_max = max(AUC_runs,[],2);
% figure; plot(trial_counts, auc_max, 'o')
%  errorbar(trial_counts, auc_mean{i}, auc_std{i}, 'ko-', 'MarkerSize', 5);
k = linspace(0,0.6, length(AUC_runs));
 color = [k(i), k(i), k(i)];

 fillOut = fill([trial_counts fliplr(trial_counts)],[(auc_mean{i}+auc_std{i})' fliplr((auc_mean{i}-auc_std{i})')],color, 'FaceAlpha', 0.1,'linestyle','none'); hold on;
 lineOut(i) = plot(trial_counts,auc_mean{i},'Color',color,'linewidth',1.5);
xlabel('number of trials'); ylabel('AUROC'); axis square; box off;
label = [num2str(cont(i)) '% contrast'];
    leg{i} = label;
hold on; 
end 
legend(lineOut, leg)
% savefig(gcf, 'AUROCvsNrTrials.fig')
%%
% Assume X is the feature matrix and y is the class labels
cv = cvpartition(y,'Holdout',0.2);
X_train = X(cv.training,:);
y_train = y(cv.training,:);
X_test = X(cv.test,:);
y_test = y(cv.test,:);

% Stratified subsampling
n_samples = 1000;
idx_pos = find(y_train == 1);
idx_neg = find(y_train == 0);
n_samples_pos = ceil(n_samples * sum(y_train == 1) / length(y_train));
n_samples_neg = ceil(n_samples * sum(y_train == 0) / length(y_train));
X_train_sub = [datasample(X_train(idx_pos,:), n_samples_pos, 'Replace', true);  datasample(X_train(idx_neg,:), n_samples_neg, 'Replace', true)];
y_train_sub = [ones(n_samples_pos,1); zeros(n_samples_neg,1)];

%%
trial_counts = [10,15,20,30,40,50,60,70,80];
n_subsamples = length(trial_counts);
Bcoeffs = zeros(5, n_subsamples);

for k = 1:5
    for i = 1:n_subsamples
        % Use datasample to randomly subsample trial_counts(i) trials with replacement
        subsampled_data{i} = datasample(tracedata{k}, trial_counts(i), 1, 'Replace', true);
        % Extract signal and noise data
         for s = 1: trial_counts(i)
            [subtrace{i}(s), peak_index_5perc_sub(s)] = getPeakPostStim(subsampled_data{i}(s,:), 15, 1, 2);
            % [subnoise{i}(s), peak_index_noise_sub(s)] = getPeakPostStim(subsampled_data{i}(s,:), 15, 0.1, 1.9);
            subnoise{i}(s)=std(subsampled_data{i}(s,24:30),0,2);
        end
        % for s = 1: trial_counts(i)
        %     % subfiveperc{i}(s) = max(subsampled_data{i}(s,31:35));
        %     subsampled_data_noise{i}(s) = max(subsampled_data{i}(s,28:30));
        % end
        signal_vec = subtrace{i}';
        noise_vec = subnoise{i}';
%         signal_data = subsampled_data(:, 1); % assuming the signal is in the first column
%         noise_data = subsampled_data(:, 2:end); % assuming the noise is in the remaining columns
        % Calculate mean and covariance matrices for signal and noise distributions
        mu_signal = mean(signal_vec);
        mu_noise = mean(noise_vec); % calculate mean across columns
        cov_signal = cov(signal_vec);
        cov_noise = cov(noise_vec);
        % Calculate Bhattacharyya coefficient
        Bcoeff = exp(-1/4 * (mu_signal - mu_noise)' * inv((cov_signal + cov_noise)/2) * (mu_signal - mu_noise)) * ...
                 sqrt(det((cov_signal + cov_noise)/2) / sqrt(det(cov_signal) * det(cov_noise)));
        % Save Bhattacharyya coefficient for this subsample size
        Bcoeffs(k, i) = Bcoeff;
    end
end

