%% load analyzed individual session's dataset and compile them 

% EK 23 

close all; clear all; clc;
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
% addpath(genpath('\\ad.gatech.edu\bme\labs\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))

fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\';
% Animal = {'M230130_2', 'M230131_1'};
Animal = {'M230220_1'};
% expID = {'26-Feb-2023_3', '27-Feb-2023'};
expID = {'24-May-2023', '25-May-2023', '26-May-2023'};
gratings = false; % false if barmapping
nrLocs = 17; % 2 for gratings (binoc/monoc)
%% compile data

[comp, nrSess] = compileAnalyzedData(fPath, Animal, expID, gratings);

%%
cont_compiled = [];
for exp = 1: size(comp,1) % find each experiment's number of contrasts
    if gratings 
        cont(exp) = length(comp(exp,1).contrasts); 
        nrContrasts = max(cont);
    else 
        cont(exp) = length(comp(exp,1).cont1);
        cont_compiled = horzcat(cont_compiled, comp(exp,1).cont1);
        cont_compiled = unique(floor(cont_compiled));
        cont_compiled = sort(cont_compiled, 'descend');
        nrContrasts = length(cont_compiled);
    end
end 

%% 
temp = cell(1,nrLocs); % exp x location
temp_gr = cell(length(expID),nrLocs); % exp x location
data = cell(size(comp,1),nrLocs); % exp x location
for exp = 1: size(comp,1) % experiments
    for s = 1: nrSess(exp) % sessions
        
        if gratings % sort binoc and monoc to cell
            %             for i =  1: length(cont)
            if cont(exp) == 2 % fill NAN for full contrast case
                comp(exp,s).avgTempTrace.binoc(cont(exp)+1: nrContrasts, :) = NaN;
                comp(exp,s).avgTempTrace.binoc([nrContrasts,cont(exp)],:) = comp(exp, s).avgTempTrace.binoc([cont(exp), nrContrasts],:);
                comp(exp,s).avgTempTrace.monoc(cont(exp)+1: nrContrasts, :) = NaN;
                comp(exp,s).avgTempTrace.monoc([nrContrasts,cont(exp)],:) = comp(exp, s).avgTempTrace.monoc([cont(exp), nrContrasts],:);
            end
            temp_gr{exp,1} = cat(3, temp_gr{exp,1}, comp(exp,s).avgTempTrace.binoc); % concatenate across sessions
            temp_gr{exp,2} = cat(3, temp_gr{exp,2}, comp(exp,s).avgTempTrace.monoc);
            
        else % barmapping
            for i = 1: length(comp(1,1).avgTempTrace) % locations
                comp(exp,s).avgTempTrace{end,i} = comp(exp,s).avgTempTrace{end,1}; % for 0 contrast
            end
            for i = 1: length(comp(1,1).avgTempTrace) % locations
                for j = 1: cont(exp)
                    temp{i} = vertcat(temp{i}, comp(exp,s).avgTempTrace{j,i});
                end
            end
            
            for i = 1: length(comp(1,1).avgTempTrace)
                temp{i}(cont(exp)+1 : nrContrasts, :) = NaN; %
                if cont(exp) == 2 % fill NAN for full contrast case
                    if comp(exp,s).cont1 == [100, 0]
                        temp{i}([nrContrasts,2],:) = temp{i}([2, nrContrasts],:);
                    else
                        idx = [];
                    for k = 1: length(comp(exp,s).cont1)
                        idx = [idx find(cont_compiled == floor(comp(exp,s).cont1(k)))];
                    end
                    temp{i}(idx,:) = temp{i}(1:cont(exp),:);
                    n = 1:nrContrasts; 
                    idx = setdiff(n, idx);                    
                    temp{i}(idx,:) = NaN;                        
                    end
                    
                elseif cont(exp) ~= 2 && cont(exp) ~= nrContrasts
                    idx = [];
                    for k = 1: length(comp(exp,s).cont1)
                        idx = [idx find(cont_compiled == floor(comp(exp,s).cont1(k)))];
                    end
                    temp{i}(idx,:) = temp{i}(1:cont(exp),:);
                    temp{i}(cont(exp),:) = NaN;
                end
            end
    
            for i = 1: length(comp(1,1).avgTempTrace) % locations
            data{exp,i} = cat(3, data{exp,i},  temp{i}); 
            end 
            temp = cell(1,nrLocs);
        end
    end
end 


if gratings
    temp = temp_gr;
else 
    temp = data; 
end 
%% mean df/f vs contrast across sessions 

avg_sessions = cellfun(@(x) nanmean(x, 3), temp, 'UniformOutput', false);
std_sessions =  cellfun(@(x) nanstd(x, 0, 3), temp, 'UniformOutput', false);

for exp = 1: size(avg_sessions,1)
    for j = 1: size(avg_sessions, 2) % locations
        for i = 1:nrContrasts
            [peak_value(i,j), peak_index] = getPeakPostStim(avg_sessions{exp,j}(i,:), 15, 1, 2);
            std_sessions_peak(i,j) = std_sessions{exp,j}(i, peak_index);
        end
    end
    % end
    
    f = figure;
    subplot(1,2,1)
    if gratings
        errorbar(comp(end,1).contrasts*100, peak_value(:,1), std_sessions_peak(:,1), 'o');
    else
        errorbar(cont_compiled, peak_value(:,5), std_sessions_peak(:,5), 'o'); % 5th bar location is binoc 
    end
    xlabel('Contrasts (%)')
    ylabel('\Deltaf/f');
    title('0 degree'); axis square; axis padded
    
    subplot(1,2,2)
    if gratings
        errorbar(comp(end,1).contrasts*100, peak_value(:,2), std_sessions_peak(:,2), 'o');
    else
        errorbar(cont_compiled, peak_value(:,12), std_sessions_peak(:,12), 'o');
    end
    xlabel('Contrasts (%)')
    ylabel('\Deltaf/f');
    title('70 degree'); axis square; axis padded
    sgtitle(['average response activity across sesssions - ' expID(exp)]);
    
    yLim = [min(min(peak_value))-max(max(std_sessions_peak)) max(max(peak_value))+max(max(std_sessions_peak))]; % Define the y-axis limit
    ax = findobj(f,'type','axes'); % Find all axes handles in the figure
    set(ax,'ylim',yLim);
    
    if gratings
        dir = [fPath filesep Animal filesep 'gratings'];
    else
        dir = [fPath filesep Animal filesep 'barmapping'];
    end
    cd([dir filesep expID{exp}])
    savefig(f, 'meandff_over_contrasts_sessions.fig')
    saveas(f, 'meandff_over_contrasts_sessions.png')

end 
%% across days 

avg_exp = cell(1,size(avg_sessions,2));
for exp =1: size(avg_sessions,1)
    for i = 1: size(avg_sessions,2) % locations 
        avg_exp{1,i} = cat(3, avg_exp{1,i}, temp{exp,i}); % concatenate over sessions
 
    end
end 

avg_days = cellfun(@(x) nanmean(x, 3), avg_exp, 'UniformOutput', false);
std_days =  cellfun(@(x) nanstd(x, 0, 3), avg_exp, 'UniformOutput', false);

for j = 1: size(avg_days, 2) % locations
    for i = 1:nrContrasts
        [peak_value_days(i,j), peak_index] = getPeakPostStim(avg_days{j}(i,:), 15, 1, 2);
        std_days_peak(i,j) = std_days{j}(i, peak_index);
    end
end 
% end 

f = figure;
subplot(1,2,1)
if gratings 
    errorbar(comp(end,1).contrasts, peak_value_days(:,1), std_days_peak(:,1), 'o');
else 
    errorbar(cont_compiled, peak_value_days(:,5), std_days_peak(:,5), 'o');
end 
xlabel('Contrasts (%)')
ylabel('\Deltaf/f');
title('Binoc Session'); axis square; axis padded

subplot(1,2,2)

if gratings 
errorbar(comp(end,1).contrasts, peak_value_days(:,2), std_days_peak(:,2), 'o');
else 
  errorbar(cont_compiled, peak_value_days(:,12), std_days_peak(:,12), 'o');
end 
xlabel('Contrasts (%)')
ylabel('\Deltaf/f');
title('Monoc Session'); axis square; axis padded
sgtitle(['average response activity across days - ' Animal]);

yLim = [min(min(peak_value_days))-max(max(std_days_peak)) max(max(peak_value_days))+max(max(std_days_peak))]; % Define the y-axis limit
ax = findobj(f,'type','axes'); % Find all axes handles in the figure
set(ax,'ylim',yLim);

if gratings
  dir = [fPath filesep Animal filesep 'gratings'];
else
  dir = [fPath filesep Animal filesep 'barmapping']; 
end
cd(dir)
savefig(f, 'meandff_over_contrasts_days.fig')
saveas(f, 'meandff_over_contrasts_days.png')
save('compiled ROI traces_days.mat', 'avg_exp', 'expID')
