close all; clear all; clc;
addpath(genpath('Y:\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))
addpath(genpath('Y:\haider\Code\behaviorAnalysis'))
% addpath(genpath('\\ad.gatech.edu\bme\labs\haider\Data\analyzedData\EKK\WFI\Codes\Analysis_EK'))


fPath = 'Y:\haider\Data\analyzedData\EKK\WFI\stats\behavior';
load([fPath filesep 'compBeh.mat']);
%%
% animal = {'M230220_1', 'M230209_1'};
animal = {'M230717_1'};

% mC = [0; 0.02; 0.05; 0.1; 0.18; 0.33];
mC = [0; 0.02; 0.03; 0.04; 0.05; 0.1; 0.18; 0.33];

% bC = [0; 0.01; 0.02; 0.05; 0.1; 0.18];
bC = [0; 0.01; 0.02; 0.03; 0.05; 0.1; 0.18; 0.33];
for i = 1:length (animal)
    for  d = 1:length(key.(animal{i}))
        for s = 1:length(key.(animal{i})(d).Sessions)
            if ~isempty(key.(animal{i})(d).Sessions(s).Out)
                comp{i}(d,s) = analyze_active_hlab_corrected_KP(key.(animal{i})(d).Sessions(s).Out, 'both');
                if ~isequal(comp{i}(d,s).mC, mC)

                    idx = find(~ismember(mC, comp{i}(d,s).mC) ==0);
                    idx1 = find(ismember(mC, comp{i}(d,s).mC) ==0);
                    temp = comp{i}(d,s).mdprime;
                    comp{i}(d,s).mdprime(end+1:length(mC)) = NaN;
                    comp{i}(d,s).mdprime(idx(2:end))= temp;
                    comp{i}(d,s).mdprime(idx1)= NaN; comp{i}(d,s).mdprime(1) = NaN;
                    clear temp
                    temp =comp{i}(d,s).mrate;
                    comp{i}(d,s).mrate(end+1:length(mC)) = NaN;
                    comp{i}(d,s).mrate(idx(1:end))= temp;
                    comp{i}(d,s).mrate(idx1)= NaN; 
                    clear temp
                    temp =comp{i}(d,s).mrt.avg;
                    comp{i}(d,s).mrt.avg(end+1:length(mC)) = NaN;
                    comp{i}(d,s).mrt.avg(idx(1:end))= temp;
                    comp{i}(d,s).mrt.avg(idx1)= NaN; 
                    temp =comp{i}(d,s).mrt.sem;
                    comp{i}(d,s).mrt.sem(end+1:length(mC)) = NaN;
                    comp{i}(d,s).mrt.sem(idx(1:end))= temp;
                    comp{i}(d,s).mrt.sem(idx1)= NaN; 
                    % comp{i}(d,s).mrt.sem = [comp{i}(d,s).mrt.sem(1) NaN comp{i}(d,s).mrt.sem(2:end)];
                    comp{i}(d,s).mC = mC;
                end
                if ~isequal(comp{i}(d,s).bC, bC)
                    idx = find(~ismember(bC, comp{i}(d,s).bC) ==0);
                    idx1 = find(ismember(bC, comp{i}(d,s).bC) ==0);
                    temp = comp{i}(d,s).bdprime;
                    comp{i}(d,s).bdprime(end+1:length(bC)) = NaN;
                    comp{i}(d,s).bdprime(idx(2:end))= temp;
                    comp{i}(d,s).bdprime(idx1)= NaN; comp{i}(d,s).bdprime(1) = NaN;
                    clear temp
                    temp =comp{i}(d,s).brate;
                    comp{i}(d,s).brate(end+1:length(bC)) = NaN;
                    comp{i}(d,s).brate(idx(1:end))= temp;
                    comp{i}(d,s).brate(idx1)= NaN; 
                    clear temp
                    temp =comp{i}(d,s).brt.avg;
                    comp{i}(d,s).brt.avg(end+1:length(bC)) = NaN;
                    comp{i}(d,s).brt.avg(idx(1:end))= temp;
                    comp{i}(d,s).brt.avg(idx1)= NaN; 
                    temp =comp{i}(d,s).brt.sem;
                    comp{i}(d,s).brt.sem(end+1:length(bC)) = NaN;
                    comp{i}(d,s).brt.sem(idx(1:end))= temp;
                    comp{i}(d,s).brt.sem(idx1)= NaN; 
                    % comp{i}(d,s).mrt.sem = [comp{i}(d,s).mrt.sem(1) NaN comp{i}(d,s).mrt.sem(2:end)];
                    comp{i}(d,s).bC = bC;
                    % comp{i}(d,s).bdprime = [ NaN comp{i}(d,s).bdprime(1:end)];
                    % comp{i}(d,s).brate = [comp{i}(d,s).brate(1) NaN comp{i}(d,s).brate(2:end)];
                    % comp{i}(d,s).brt.avg = [comp{i}(d,s).brt.avg(1) NaN comp{i}(d,s).brt.avg(2:end)];
                    % comp{i}(d,s).brt.sem = [comp{i}(d,s).brt.sem(1) NaN comp{i}(d,s).brt.sem(2:end)];
                    % comp{i}(d,s).bC = bC;
                end
            end
        end
    end
end

%%
mdprime = []; bdprime =[];
mrate = []; brate = [];
brt.avg = []; brt. sem = [];
mrt.avg = []; mrt. sem =[];

for i = 1:length (animal)
    for  d = 1:length(key.(animal{i}))
        for s = 1:length(key.(animal{i})(d).Sessions)
            if ~isempty(comp{i}(d,s).mC)
                mdprime = [mdprime; comp{i}(d,s).mdprime];
                bdprime = [bdprime; comp{i}(d,s).bdprime];
                brt.avg = [brt.avg; comp{i}(d,s).brt.avg];
                mrt.avg = [mrt.avg; comp{i}(d,s).mrt.avg];
                brt.sem = [brt.sem; comp{i}(d,s).brt.sem];
                mrt.sem = [mrt.sem; comp{i}(d,s).mrt.sem];
                mrate = [mrate; comp{i}(d,s).mrate];
                brate = [brate; comp{i}(d,s).brate];
            end
        end
    end
end

activeOut.mC = mC; activeOut.bC = bC;
activeOut.mdprime = nanmean(mdprime,1);
activeOut.bdprime = nanmean(bdprime,1);
activeOut.brt.avg = nanmean(brt.avg,1);
activeOut.mrt.avg = nanmean(mrt.avg,1);
activeOut.brt.sem = nanmean(brt.sem,1);
activeOut.mrt.sem = nanmean(mrt.sem,1);
activeOut.mrate = nanmean(mrate,1);
activeOut.brate = nanmean(brate,1);

%%
SPs = [5, 5, 5, 5; 0.01, 0.05, 5, 1; -5, -5, -5, -5];   %
figure; 
subplot 321 % m dprime
errorbar (activeOut.mC(2:end), activeOut.mdprime(2:end), nanstd(mdprime(:,2:end),[],1)/sqrt(size(mdprime(:,2:end),1)), '-', 'CapSize',0, 'LineWidth', 1.5); hold on; 
hold on; 
if length(activeOut.mC(2:end)) >= 4
    [ffit_monoc, curve_monoc] = FitPsycheCurveWH([activeOut.mC(2:end)], activeOut.mdprime(2:end), SPs);
    plot(curve_monoc(:,1), curve_monoc(:,2), 'Color', [1,0,0,0.75],'LineWidth', 1); 
end
setfig
title(['Monocular'])
ylabel(['d prime'])
yline([0], '--k')
xlim([-0.05 0.5])
ylim([-0.5 3])
axis square
xlabel(['Contrast'])


if isempty(activeOut.mC)
    c50_m = NaN;
else

subplot 323
SPs = [5, 5, 5, 5; 0.01, 0.05, 5, 1; -5, -5, -5, -5];
errorbar (activeOut.mC(2:end), activeOut.mrate(2:end), nanstd(mrate(:,2:end),[],1)/sqrt(size(mrate,1)), '-', 'CapSize',0,'LineWidth', 1.5); hold on; 
 hold on; 
 errorbar(activeOut.mC(1), activeOut.mrate(1), nanstd(mrate(:,1),[],1)/sqrt(size(mrate,1)), 'ko', 'CapSize',0,'LineWidth', 1.5); hold on; 
if length(activeOut.mC(2:end)) >= 4
    [ffit_monoc, curve_monoc] = FitPsycheCurveWH([activeOut.mC(1:end)], activeOut.mrate(1:end), SPs);
    plot(curve_monoc(:,1), curve_monoc(:,2), 'Color', [1,0,0,0.75], 'LineWidth', 1); 
end


far = activeOut.mrate(1);
lapse = activeOut.mrate(end);

if length(activeOut.mC(2:end)) >= 4
hh = (lapse - far)./2
p_c50_m = hh+far;
[~,y] = min(abs(curve_monoc(:,2)-p_c50_m));

c50_m = curve_monoc(y,1);
else c50_m = 0; p_c50_m = 0;
end

setfig
ylim([-.05 1])
xlim([-0.05 0.5])
ylabel(['hit rate'])
axis square
xlabel(['Contrast'])
if isnan(c50_m)
    title(['too easy: all > 50% hit rate'])
    yline([0.5], '--k');
else
title(['threshold = ', num2str(round(c50_m,2))])
yline([p_c50_m], '--k');
xline([c50_m], '--k');
end


subplot 325
errorbar(activeOut.mC(2:end), activeOut.mrt.avg(2:end), activeOut.mrt.sem(2:end), '-ok', 'CapSize', 0, 'LineWidth', 1.5); hold on; 
setfig
xlim([-0.05 0.5])
ylabel(['reaction time (s)'])
axis square
xlabel(['Contrast'])
end

if isempty(activeOut.bC)
    c50_b = NaN;
else
subplot 322
errorbar (activeOut.bC(2:end), activeOut.bdprime(2:end), nanstd(bdprime(:,2:end),[],1)/sqrt(size(bdprime(:,2:end),1)), '-', 'CapSize',0, 'LineWidth', 1.5); hold on; 

if length(activeOut.bC(2:end)) >= 4
    [ffit_binoc, curve_binoc] = FitPsycheCurveWH([activeOut.bC(2:end)], activeOut.bdprime(2:end), SPs);
    plot(curve_binoc(:,1), curve_binoc(:,2), 'Color', [1,0,0,0.75], 'LineWidth', 1); 
end
setfig
title(['Binocular'])
ylabel(['d prime'])
yline([0], '--k')
xlim([-0.05 0.5])
ylim([-0.5 3])
axis square
xlabel(['Contrast'])

subplot 324
SPs = [5, 5, 5, 5; 0.01, 0.05, 5, 1; -5, -5, -5, -5];
errorbar (activeOut.bC(2:end), activeOut.brate(2:end), nanstd(brate(:,2:end),[],1)/sqrt(size(brate,1)), '-', 'CapSize',0,'LineWidth', 1.5); hold on; 
 hold on; 
 errorbar(activeOut.bC(1), activeOut.brate(1), nanstd(brate(:,1),[],1)/sqrt(size(brate,1)), 'ko', 'CapSize',0,'LineWidth', 1.5); hold on; 
if length(activeOut.bC(2:end)) >= 4
    [ffit_binoc, curve_binoc] = FitPsycheCurveWH([activeOut.bC(2:end)], activeOut.brate(2:end), SPs);
    plot(curve_binoc(:,1), curve_binoc(:,2), 'Color', [1,0,0,0.75], 'LineWidth', 1); 
end

far = activeOut.brate(1);
lapse = activeOut.brate(end);

if length(activeOut.bC(2:end)) >= 4
    hh = (lapse - far)./2;
    p_c50_b = hh+far;
    [~,y] = min(abs(curve_binoc(:,2)-p_c50_b));

    c50_b = curve_binoc(y,1);
else c50_b = 0; p_c50_b = 0;
end


setfig
ylim([-.05 1])
xlim([-0.05 0.5])
ylabel(['hit rate'])
axis square
xlabel(['Contrast'])
if isnan(c50_b)
    title(['too easy: all > 50% hit rate'])
    yline([0.5], '--k');
else
title(['threshold = ', num2str(round(c50_b,2))])
yline([p_c50_b], '--k');
xline([c50_b], '--k');
end

subplot 326
errorbar(activeOut.bC(2:end), activeOut.brt.avg(2:end), activeOut.brt.sem(2:end), '-ok', 'CapSize', 0, 'LineWidth', 1.5); hold on; 
setfig
xlim([-0.05 0.5])
ylabel(['reaction time (s)'])
axis square
xlabel(['Contrast'])
end