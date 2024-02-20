
cd('Y:\haider\Data\analyzedData\EKK\WFI\stats\gratings')
passive = load ('compiled ROI traces_animals.mat');
cd('Y:\haider\Data\analyzedData\EKK\WFI\stats\behavior')
active = load('compiled ROI traces_animals_hit.mat');

avg_animal_p = cellfun(@(x) nanmean(x, 3), passive.temp, 'UniformOutput', false);
std_animal_p =  cellfun(@(x) nanstd(x, 0, 3), passive.temp, 'UniformOutput', false);
nrContrasts = size(passive.temp{1},1);

for j = 1: size(avg_animal_p, 2) % locations
    for i = 1:nrContrasts
        [peak_value_p(i,j), peak_index] = getPeakPostStim(avg_animal_p{j}(i,:), 15, 1, 2,0);
        std_peak_p(i,j) = std_animal_p{j}(i, peak_index);
    end
end 

avg_animal_a = cellfun(@(x) nanmean(x, 3), active.temp, 'UniformOutput', false);
std_animal_a =  cellfun(@(x) nanstd(x, 0, 3), active.temp, 'UniformOutput', false);
nrContrasts_a = length(active.data(1).idxb);
for j = 1: size(avg_animal_a, 2) % locations
    for i = 1:nrContrasts_a
        [peak_value_a(i,j), peak_index] = getPeakPostStim(avg_animal_a{j}(i,:), 15, 1, 2,0);
        std_peak_a(i,j) = std_animal_a{j}(i, peak_index);
    end
end 

%% 
m_p =peak_value_p(1:5,2); % 0 2 5 10 33;
m_a = peak_value_a([1:4, 6], 2);
m_p_std = std_peak_p(1:5,2);
m_a_std = std_peak_a([1:4, 6], 2);

% Create x-axis labels for each group
x_labels = {'passive', 'behavior'};

% Combine the mean values and standard deviations for both groups
all_means = [m_p; m_a];
all_stds = [m_p_std; m_a_std];

% Create x-axis values
x = 1:length(x_labels);

% Plot the data points for each contrast and group
% figure;
% hold on;
% scatter(repmat(x(1), 1, length(m_p)), m_p, 'ro');
% scatter(repmat(x(2), 1, length(m_a)), m_a, 'bo');
% hold off;

% Add error bars for each group
figure; 
k = linspace(0.7,0.1, 5);
for i = 2:5
    plot([1 2], [m_p(i)*100 m_a(i)*100], 'LineWidth', 1.5,'color', [0.8 0.8 0.8]); hold on
end
for i = 2:5
t(i) = errorbar(1, m_p(i)*100, m_p_std(i)*100, 'o', 'LineWidth', 1, 'CapSize',0, 'color', [k(i), k(i), k(i)],'MarkerFaceColor',[k(i), k(i), k(i)]);
errorbar(2, m_a(i)*100, m_a_std(i)*100, 'o', 'LineWidth', 1, 'CapSize',0, 'color', [0.77,0.44,0.44],'MarkerFaceColor',[0.77,0.44,0.44]);
hold on
end 
% errorbar([1 2], [m_p(i) m_a(i)], [m_p_std(i) m_a_std(i)], 'o', 'LineWidth', 1.5, 'CapSize',0, 'color', [0.50,0.50,0.50], 'MarkerFaceColor',[0.20,0.20,0.20]);

% Add labels and legend

set(gca, 'XTick', x, 'XTickLabel', x_labels);
axis padded; axis square; box off
set(findobj(gcf, 'type', 'axes'), 'xlim', [0 3])
legend(t(2:end) , '2%', '5%', '10%', '33%'); legend boxoff
ylabel('Mean Peak Value (%)');

%%
for ii = 1:2
for j = 1: size(passive.temp{1,ii}, 3) % locations
    for i = 1:nrContrasts
        [a{ii}(i,j), peak_index] = getPeakPostStim(passive.temp{ii}(i,:,j), 15, 1, 2,0);
        
        % std_peak_p(i,j) = std_animal_p{j}(i, peak_index);
    end
end 
end

for ii = 1:2
for j = 1: size(passive.temp{1,ii}, 3) % locations
    for i = 1:6
   [aa{ii}(i,j), peak_index] = getPeakPostStim(active.temp{ii}(i,:,j), 15, 1, 2,0);
        
        % std_peak_p(i,j) = std_animal_p{j}(i, peak_index);
    end
end 
end

% Perform paired t-test between the groups for each contrast level
p_values = zeros(1, length(m_p));
for i = 2:length(m_p)
    [~, p_values(i)] = ttest(single(aa{2}(i,:))', a{2}(i,:)');
end

% Add p-values as text annotations
p_value_threshold = 0.05;
for i = 2:5
    if ~isnan(p_values(i)) && p_values(i) < p_value_threshold
        text(x(2) + 0.1, m_a(i,1)*100, "*", 'FontSize', 12); hold on;
    end
end
text(x(1), m_p + 0.005, num2str(p_values', '%0.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'r');
text(x(2), m_a + 0.005, num2str(p_values', '%0.4f'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'b');
%%

percent_increase = 100*((m_a - m_p)./m_p); 
figure;
 bar ([2 5 10 33], percent_increase(2:end))
axis padded; axis square; box off
set(findobj(gcf, 'type', 'axes'), 'xlim', [-5 40])
xlabel('Contrast (%)', 'FontSize', 12); ylabel('Percent Increase (%)', 'FontSize', 12)
 %%
f = figure;
subplot(1,2,1)
 
errorbar(analyzed_comp(end,1).contrasts*100, peak_value(:,1)*100, std_peak(:,1)*100, 'ko' ,'LineWidth', 1.5);
xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('Binocular V1'); axis square; axis padded; box off

subplot(1,2,2)
if gratings 
errorbar(analyzed_comp(end,1).contrasts*100, peak_value(:,2)*100, std_peak(:,2)*100, 'ko', 'LineWidth', 1.5);
else 
    errorbar(cont, peak_value(:,12), std_peak(:,12), 'o');
end
xlabel('Contrasts (%)')
ylabel('\Deltaf/f (%)');
title('Monocular V1'); axis square; axis padded; box off
sgtitle(['average response activity across animals (N = ' num2str(length(Animal)) ')']);
