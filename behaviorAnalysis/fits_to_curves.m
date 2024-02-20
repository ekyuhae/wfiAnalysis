function [c50_m, c50_b] = fits_to_curves(activeOut)

SPs = [5, 5, 5, 5; 0.01, 0.05, 5, 1; -5, -5, -5, -5];   % 
figure; 
subplot 321
plot(activeOut.mC(2:end), activeOut.mdprime, '-o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'); hold on; 
if length(activeOut.mC(2:end)) >= 4
    [ffit_monoc, curve_monoc] = FitPsycheCurveWH([activeOut.mC(2:end)], activeOut.mdprime, SPs);
    plot(curve_monoc(:,1), curve_monoc(:,2), 'Color', [1,0,0,0.75],'LineWidth', 2); 
end
setfig
title(['Monocular'])
ylabel(['d prime'])
yline([0], '--k')
xlim([-0.05 0.85])
ylim([-0.5 3])
axis square
xlabel(['Contrast'])
% 
% subplot 323
% SPs = [5, 5, 5, 5; 0.01, 0.05, 5, 1; -5, -5, -5, -5];
% [ffit_monoc, curve_monoc] = FitPsycheCurveWH([activeOut.mC(2:end)], activeOut.mrate(2:end), SPs);
% plot(activeOut.mC(2:end), activeOut.mrate(2:end), '-o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'); hold on; 
% plot(activeOut.mC(1), activeOut.mrate(1), 'o', 'MarkerFaceColor', [128 128 128]./255, 'MarkerEdgeColor', 'k'); hold on; 
% plot(curve_monoc(:,1), curve_monoc(:,2), 'Color', [1,0,0,0.75], 'LineWidth', 2); 
% [~,y] = min(abs(curve_monoc(:,2)-0.5));
% if isempty(y)
%    c50_m = NaN;
% else
% c50_m = curve_monoc(y,1);
% end
% setfig
% ylim([-.05 1])
% xlim([-0.05 0.85])
% ylabel(['hit rate'])
% axis square
% xlabel(['Contrast'])
% if isnan(c50_m)
%     title(['too easy: all > 50% hit rate'])
%     yline([0.5], '--k');
% else
% title(['C_5_0 = ', num2str(c50_m)])
% yline([0.5], '--k');
% xline([c50_m], '--k');
% end

%% testing
if isempty(activeOut.mC)
    c50_m = NaN;
else

subplot 323
SPs = [5, 5, 5, 5; 0.01, 0.05, 5, 1; -5, -5, -5, -5];
plot(activeOut.mC(2:end), activeOut.mrate(2:end), '-o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'); hold on; 
plot(activeOut.mC(1), activeOut.mrate(1), 'o', 'MarkerFaceColor', [128 128 128]./255, 'MarkerEdgeColor', 'k'); hold on; 
if length(activeOut.mC(2:end)) >= 4
    [ffit_monoc, curve_monoc] = FitPsycheCurveWH([activeOut.mC(1:end)], activeOut.mrate(1:end), SPs);
    plot(curve_monoc(:,1), curve_monoc(:,2), 'Color', [1,0,0,0.75], 'LineWidth', 2); 
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
xlim([-0.05 0.85])
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
xlim([-0.05 0.85])
ylabel(['reaction time (s)'])
axis square
xlabel(['Contrast'])
end

if isempty(activeOut.bC)
    c50_b = NaN;
else
subplot 322
plot(activeOut.bC(2:end), activeOut.bdprime, '-o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'); hold on; 
if length(activeOut.bC(2:end)) >= 4
    [ffit_binoc, curve_binoc] = FitPsycheCurveWH([activeOut.bC(2:end)], activeOut.bdprime, SPs);
    plot(curve_binoc(:,1), curve_binoc(:,2), 'Color', [1,0,0,0.75], 'LineWidth', 2); 
end
setfig
title(['Binocular'])
ylabel(['d prime'])
yline([0], '--k')
xlim([-0.05 0.85])
ylim([-0.5 3])
axis square
xlabel(['Contrast'])

subplot 324
SPs = [5, 5, 5, 5; 0.01, 0.05, 5, 1; -5, -5, -5, -5];
plot(activeOut.bC(2:end), activeOut.brate(2:end), '-o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k'); hold on; 
plot(activeOut.bC(1), activeOut.brate(1), 'o', 'MarkerFaceColor', [128 128 128]./255, 'MarkerEdgeColor', 'k'); hold on; 
if length(activeOut.bC(2:end)) >= 4
    [ffit_binoc, curve_binoc] = FitPsycheCurveWH([activeOut.bC(2:end)], activeOut.brate(2:end), SPs);
    plot(curve_binoc(:,1), curve_binoc(:,2), 'Color', [1,0,0,0.75], 'LineWidth', 2); 
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
xlim([-0.05 0.85])
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
xlim([-0.05 0.85])
ylabel(['reaction time (s)'])
axis square
xlabel(['Contrast'])
end