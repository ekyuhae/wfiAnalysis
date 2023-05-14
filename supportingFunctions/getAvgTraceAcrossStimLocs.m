function [mean, std] = getAvgTraceAcrossStimLocs(avg_animal, contrasts, fr)
% each column of avg_animal cell array (or avgTempTrace) indicates stimulus location 
% fr = original frame rate 
% contrasts = [100 ... 0]
% EK 23 


temp = [];
for i = 1:length(avg_animal)
    temp = cat(3, temp, avg_animal{i});
end 

avg_locs = nanmean(temp, 3);

fr = fr/2;
t = 0:1/fr:(length(avg_animal{1,1})-1)/fr; % convert frame into timescale 
t = t - t(fr*2+1); 
% t = t - t(floor(fr*0.5)+1);
% k = linspace(0.9,0, length(contrasts));
k = linspace(0,0.9, length(contrasts));

% figure; 
% % stdshade(avg_locs)
% for i = 1: size(avg_locs,1)
% 
%     legend_names{i} = sprintf('%d%%', contrasts(i));
%     color = [k(i), k(i), k(i)];
%     plot(t, avg_locs(i,:)*100, 'Color', color);
%     xlim([t(15) t(end)])
%     hold on
% box off; axis square
% 
% % stdshade(temp1(:,:,i)); hold on;
% end 
%   title('avg response across stimulus locations');
%         xline(t(fr*2+1), ':k'); hold off  
%    legend(legend_names); legend('boxoff');
%    xlabel ('relative time to stimulus onset (s)')
%    ylabel('dF/F (%)')
%    


%%
figure;

for i = 1: size(avg_locs,1)
    
    legend_names{i} = sprintf('%.2f%%', contrasts(i));
    color = [k(i), k(i), k(i)];
    [mean(i,:), std(i,:), lineout(i), shadeout(i)] = stdshade(temp(i,:,:), 0.1, color, t, [] ); hold on;
    xlim([t(15) t(end)]) %0410
    hold on
   
end
title('avg response across stimulus locations'); box off; axis square
xline(t(fr*2+1), ':k'); hold off
legend(lineout, legend_names); legend('boxoff');
xlabel ('relative time to stimulus onset (s)')
ylabel('dF/F')

% savefig(gcf, 'meandff_over_locations_animals.fig')
% saveas(gcf, 'meandff_over_locations_animals.png')