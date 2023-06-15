% perform linear discriminant analysis on df/f data during barmapping
% experiment 

close all; clear all; clc;
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\';
Animal = {'M221118_1','M230130_2','M230131_1'};
cont = [100 75 50 25 5 0];
full = [1 0 0 0 0 1]; five = [0 0 0 0 1 1]; multi = [1 1 1 1 0 1]; % 1= true for each contrast value

data = struct();
data.(Animal{1})(1,:) = {'02-Feb-2023_1', '03-Feb-2023', '08-Feb-2023_2', '15-Feb-2023_2', '18-Feb-2023', '09-Mar-2023'};
% data.(Animal{1})(1,:) = {'15-Feb-2023_2', '18-Feb-2023', '09-Mar-2023'};
data.(Animal{1})(2,:) ={full, full, full, multi, multi, five};
data.(Animal{2})(1,:) = {'20-Feb-2023', '24-Feb-2023', '09-Mar-2023'};
data.(Animal{2})(2,:) ={full, multi, five};
data.(Animal{3})(1,:) = {'01-Mar-2023', '01-Mar-2023_1'};
data.(Animal{3})(2,:) ={full, multi};

combinedROIs = true; % picking [best-1 best best+1] ROI response

%% data extraction and formatting 
traces = struct();
for i = 1:length(Animal)
   for exp= 1:numel(data.(Animal{i})(1,:))
    temp  = compileIndividTrial(Animal{i}, data.(Animal{i}){1,exp}, fPath, sum(data.(Animal{i}){2,exp}));
    % dFFtraces = [dFFtraces; temp];
    traces(i,exp).trace = temp;
   end
end

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
%% data transformation
if combinedROIs
    trace5 = [];
    trace25 = [];
    trace50 = [];
    trace75 = [];
    trace100 = [];

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

else % original ROI traces

    trace5 = [];
    trace25 = [];
    trace50 = [];
    trace75 = [];
    trace100 = [];

    for i = 1:17
        trace5 = [trace5; dFFtraces{5,i}];
        trace25 = [trace25; dFFtraces{4,i}];
        trace50 = [trace50; dFFtraces{3,i}];
        trace75 = [trace75; dFFtraces{2,i}];
        trace100 = [trace100; dFFtraces{1,i}];
    end
    tracedata{5} = trace5; tracedata{4} = trace25; tracedata{3} = trace50; tracedata{2} = trace75; tracedata{1} = trace100;
    noisetrace = dFFtraces{6,1};

end
clear trace5 trace25 trace50 trace75 trace100
%% LDA for whole trials
for k = 1:length(tracedata)
    for i = 1: size(tracedata{k},1)
        tracedata1{k}(i,:) = tracedata{k}(i,:) - min(tracedata{k}(i,:));
        [signal{k}(i), peak_index_5perc(i)] = getPeakPostStim(tracedata1{k}(i,:), 15,1, 2);
        noise{k}(i) = getPeakPostStim(tracedata1{k}(i,:), 15,0.5, 1.5);
    end

    % figure; histogram(signal{k},50,'Normalization', 'probability'); hold on; histogram(noise_fromstd{k}',50,'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha',0.5)
    % legend([num2str(cont(k)) '% contrast'] , '0% contrast'); box off; legend('boxoff'); axis square
    % title ([num2str(cont(k)) '% contrast (' num2str(length(signal{k})) ' trials)']); xlabel('\DeltaF/F_0'); ylabel ('probability')

    % leg= cell(5,1);

    signal_vec = signal{k}';
    noise_vec = noise{k}';

    X = [signal_vec; noise_vec]; % predictor variable
    Y = logical([repmat(1,length(signal_vec),1); repmat(0,length(noise_vec), 1)]); %index indicating signal (1) and noise (0) -response
    n = length(Y);
    fold = 10;
    c = cvpartition(n, 'KFold', fold);

    for i = 1:fold
        X_train = X(c.training(i),:); % training set w observed data or input
        Y_train = Y(c.training(i),:); % training set w response
        X_test = X(c.test(i),:);
        Y_test = Y(c.test(i),:);

        LDA = fitcdiscr(X_train,Y_train); %apply LDA on training set
        Y_pred = predict(LDA, X_test); % apply LDA on testing set
        acc_whole(i) = sum(Y_pred == Y_test)/length(Y_test);
        precision(i) = sum(Y_pred(Y_test==1)==1) / sum(Y_pred==1);
        recall(i) = sum(Y_pred(Y_test==1)==1) / sum(Y_test==1);
        f1score(i) = 2 * (precision(i)  * recall(i)) / (precision(i)  + recall(i));
    end

    meanlda.acc(k) = mean(acc_whole); sdlda.acc(k) = std(acc_whole);
    meanlda.pre(k) = mean(precision); sdlda.pre(k) = std(precision);
    meanlda.rec(k) = mean(recall); sdlda.rec(k) = std(recall);
    meanlda.f1(k) = mean(f1score); sdlda.f1(k) = std(f1score);
    clear signal_vec noise_vec X Y c LDA Y_pred
end

figure;
subplot 221
errorbar(cont(1:5), meanlda.acc, sdlda.acc, 'ko:', LineWidth=2); xlabel('Contrast Level (%)'); ylabel('Discrimination Accuracy(%)')
% axis square; box off;
subplot 222
errorbar(cont(1:5), meanlda.pre, sdlda.pre, 'ko:', LineWidth=2); xlabel('Contrast Level (%)'); ylabel('Precision(%)')
subplot 223
errorbar(cont(1:5), meanlda.rec, sdlda.rec, 'ko:', LineWidth=2); xlabel('Contrast Level (%)'); ylabel('recall(%)')
subplot 224
errorbar(cont(1:5), meanlda.f1, sdlda.f1, 'ko:', LineWidth=2); xlabel('Contrast Level (%)'); ylabel('f1score(%)')
f = findobj(gcf, 'type','axes'); axis(f, 'square'); box(f, 'off')
% title(f,"LDA performance vs. Contrasts")
clear precision recall f1score mean_acc X Y X_train X_test Y_train Y_test Y_pred n

%% LDA on subsampled trials 
trial_counts = [15,20,25,30,40,50,60,70,85,100];
n_subsamples = length(trial_counts);
iter = [10, 15, 20, 30];
K= 5; 

for k = 1:length(tracedata)
    %%
    for i = 1:n_subsamples
        % subsampled_data = cell(n_subsamples, 1);
        % if k == 5
        for c_iteration = 1:length(iter)
        for run = 1: iter(c_iteration)

            subsampled_data = datasample(tracedata1{k}, trial_counts(i),1, 'Replace', true);

            for s = 1: trial_counts(i)
                [subtrace(s), peak_index_5perc_sub(s)] = getPeakPostStim(subsampled_data(s,:), 15, 1, 2);
                subnoise(s)= getPeakPostStim(subsampled_data(s,:), 15, 0.5, 1.5);
                % subnoise(s)=std(subsampled_data(s,24:30),0,2);
            end

            signal_vec = subtrace';
            noise_vec = subnoise';
            X = [signal_vec; noise_vec]; % predictor variable
            Y = logical([repmat(1,length(signal_vec),1); repmat(0,length(noise_vec), 1)]); %index indicating signal (1) and noise (0)
            n = length(Y);
            c = cvpartition(n,'KFold',K);
            for j = 1: K % k-fold cross validation 
                X_train = X(c.training(j),:); % predictor
                Y_train = Y(c.training(j));
                X_test = X(c.test(j),:);
                Y_test = Y(c.test(j));

                lda = fitcdiscr(X_train,Y_train); %apply LDA on training set
                Y_pred = predict(lda, X_test);

                acc_kfold(j)= sum(Y_pred == Y_test)/length(Y_test);
                % Save accuracy for this subsample size

            end
            avgKacc(run)= mean(acc_kfold); % mean accruacy over kfolds in one subsample set
        end
        acc_sub{k}(i,c_iteration) = mean(avgKacc);
        clear avgKacc X Y X_train X_test Y_train Y_test Y_pred n lda subtrace subnoise
        end
    end
     
end 

%%
for k = 1: length(acc_sub)
    sub(k,:) = mean(acc_sub{1,k}, 2); 
    sub_sd(k,:)  = std(acc_sub{1,k}, 0,2); 
end 


%%
figure;
for i = 1:length(acc_sub)
    k = linspace(0,0.6, 5);
    color = [k(i) k(i) k(i)];
    % subplot(121)
    % errorbar(trial_counts,  mean_accs(2,:), sd_accs(2,:), ':o', 'Color', [k(1) k(1) k(1)], 'LineWidth',1.5);
    % subtitle('100% contrast'); ylabel('Discrimination Acccuracy'); xlabel('Number of Trials');
    % subplot(122)
    errorbar(trial_counts, sub(i,:), sub_sd(i,:), ':o', 'Color',[k(i) k(i) k(i)], 'LineWidth',1.5); hold on
    
    % fillOut = fill([trial_counts fliplr(trial_counts)],[(sub(i,:)+sub_sd(i,:))'; fliplr((sub(i,:)-sub_sd(i,:))')]',color, 'FaceAlpha', 0.1,'linestyle','none'); hold on;
    % lineOut(i) = plot(trial_counts,sub(i,:),'Color',color,'linewidth',1.5);
    ylabel('Discrimination Acccuracy'); xlabel('Number of Trials');
    f = findobj(gcf, 'type','axes'); axis(f, 'square'); box(f, 'off')
    % subtitle('5% contrast')
    legend({ '100%', '75%', '50%', '25%', '5%'})
end
% 
% figure; 
%    errorbar(trial_counts,  sub(5,:), sub_sd(2,:), ':o', 'Color', [k(1) k(1) k(1)], 'LineWidth',1.5);
%      ylabel('Discrimination Acccuracy'); xlabel('Number of Trials');
%     f = findobj(gcf, 'type','axes'); axis(f, 'square'); box(f, 'off')
%     ylim(f,[0.45 0.65])
%     criter = median(sub(5,:));
%%
% Perform LDA on the data
class = classify(X,Y,'linear');

% Plot the data points
figure;
gscatter(signal_vec,noise_vec,ones(size(signal_vec,1),1),'rb','o+');
xlabel('Feature 1');
ylabel('Feature 2');
title('LDA Classification');

% Plot the decision boundary
hold on
K = lda.Coeffs(1,2).Const;
L = lda.Coeffs(1,2).Linear;
f = @(x1,x2) K + [x1,x2]*L';
h = fimplicit(f,[min(X(:,1))-1,max(X(:,1))+1,min(X(:,2))-1,max(X(:,2))+1]);
set(h,'Color','k','LineWidth',2)
legend('Noise','Signal','Decision Boundary')
hold off

figure

scatter(score(:,1), score(:,2), resp)

% Compute the LDA classification boundaries
k = size(LDA.Coeffs, 2); % number of classes
if k == 2 % for binary classification
    % Compute the linear discriminant function coefficients
    coef = LDA.Coeffs(1, 2).Linear;
    % Compute the decision boundary
    b = LDA.Coeffs(1, 2).Const;
    x = [min(score(:,1)) max(score(:,1))];
    y = (-1/coef(2))*(coef(1)*x + b);
    % Plot the decision boundary
    hold on
    plot(x,y,'k','LineWidth',2)
    legend('Noise','Signal','Decision boundary','Location','Best')
    hold off
end
%% LDA on subsampled trials 
trial_counts = [15,20,30,40,55,75,100];
n_subsamples = length(trial_counts);

for k = 1:length(tracedata)
    for i = 1:n_subsamples
        % subsampled_data = cell(n_subsamples, 1);
        % if k == 5
            for run = 1:10000
                % Loop over the number of subsamples
                % Use datasample to randomly subsamplene trial_counts(i) trials with replacement
                subsampled_data = datasample(tracedata1{k}, trial_counts(i),1, 'Replace', true);
                % subsampled_data_noise{i} = subsampled_data{i}(:,15:29);
                % subsampled_data_noise_mean{i} = mean(subsampled_data_noise{i},1);
                % subsampled_data_noise{i} = getPeakPostStim(cattrace(i,:), 15, 0.5, 1.5);
                for s = 1: trial_counts(i)
                    [subtrace(s), peak_index_5perc_sub(s)] = getPeakPostStim(subsampled_data(s,:), 15, 1, 2);
                    subnoise(s)= getPeakPostStim(subsampled_data(s,:), 15, 0.5, 1.5);
                    % subnoise(s)=std(subsampled_data(s,24:30),0,2);
                end
                signal_vec = subtrace';
                noise_vec = subnoise';
                X = [signal_vec; noise_vec]; % predictor variable
                Y = logical([repmat(1,length(signal_vec),1); repmat(0,length(noise_vec), 1)]); %index indicating signal (1) and noise (0)
                n = length(Y);
                c = cvpartition(n,'HoldOut',0.2);
                X_train = X(c.training,:); % predictor
                Y_train = Y(c.training);
                X_test = X(c.test,:);
                Y_test = Y(c.test);

                lda = fitcdiscr(X_train,Y_train); %apply LDA on training set
                Y_pred = predict(lda, X_test);

                acc= sum(Y_pred == Y_test)/length(Y_test);
                % Save accuracy for this subsample size
                accs_sub{k}(i, run) = acc;
            end
             mean_accs(k,:) = meanlda(accs_sub{k},2);
        sd_accs(k,:) = std(accs_sub{k},0,2);
    end
end