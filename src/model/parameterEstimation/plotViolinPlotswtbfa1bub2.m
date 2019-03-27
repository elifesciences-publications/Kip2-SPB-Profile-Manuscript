clear
close all

myFolder = fileparts(mfilename('fullpath'));

exportPlots = true;

parameterEstimateResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'simulationResults' filesep 'parameterEstimates'];
parameterEstimateFigureResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'figures' filesep 'model'];

wtData = load([parameterEstimateResultFolder filesep 'wt_params.mat']);
bfa1bub2Data = load([parameterEstimateResultFolder filesep 'bfa1bub2_params.mat']);
prefix = 'comparison-wt-bfa1bub2';

pMin = min(wtData.parameterCombinations');
pMax = max(wtData.parameterCombinations');

%% Plot beeswarm plots
figure;
currAx = gca();
colors = currAx.ColorOrder;

concentrations.wt   = (wtData.plotSamples(1, :) ./ 140 .* 61);
concentrations.bfa1bub2 = (bfa1bub2Data.plotSamples(1, :) ./ 140 .* 61);
[xConc, gConc] = addBoxPlotGroup([], [], concentrations.wt);
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, concentrations.bfa1bub2);

violins = violinplot(xConc, gConc, 'Bandwidth', 30, 'Width', 0.4, 'ShowData', false);
violins(2).ViolinColor = colors(2, :);
xlim([0.5 2.5]);
plot([0.6 2.4], [1 1] .* pMin(1) ./ 140 .* 61, '--k');
plot([0.6 2.4], [1 1] .* pMax(1) ./ 140 .* 61, '--k');
ylabel('[Kip2]_{total} (nM)');
ylim([0 230]);

currAx.XTickLabel = {'WT', '\it{bfa1{\Delta}bub2{\Delta}}'};
currAx.XLabel.Color = [0 0 0];
currAx.YLabel.Color = [0 0 0];
currAx.XColor = [0 0 0];
currAx.YColor = [0 0 0];
xtickangle(45)
applyPaperFormatting

currFig = gcf();
currFig.Position(3) = 0.4*currFig.Position(3);
currFig.Color = 'white';

fprintf('Mean [Kip2]: %g,\tStd:%g\n', mean(concentrations.wt), std(concentrations.wt));
fprintf('Median [Kip2]: %g\n', median(concentrations.wt));

fprintf('Mean [Kip2 bfa1bub2]: %g,\tStd:%g\n', mean(concentrations.bfa1bub2), std(concentrations.bfa1bub2));
fprintf('Median [Kip2 bfa1bub2]: %g\n', median(concentrations.bfa1bub2));

if exportPlots
    print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_conc.pdf'], '-dpdf');
    export_fig([parameterEstimateFigureResultFolder filesep 'png' filesep prefix '_violin_conc.png'], '-q101', '-r300');
end
%%

k_on.wt   = (wtData.plotSamples(2, :));
k_on.bfa1bub2 = (bfa1bub2Data.plotSamples(2, :));
figure;

xConc = [];
gConc = [];
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(k_on.wt));
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(k_on.bfa1bub2));

fprintf('Mean k_on: %g,\t%g\n', mean(k_on.wt), std(k_on.wt));
fprintf('Median k_on: %g\n', median(k_on.wt));

fprintf('Mean k_on bfa1bub2del: %g\t%g\n', mean(k_on.bfa1bub2), std(k_on.bfa1bub2));
fprintf('Median k_on bfa1bub2del: %g\n', median(k_on.bfa1bub2));



k_in.wt   = (wtData.plotSamples(3, :));
k_in.bfa1bub2 = (bfa1bub2Data.plotSamples(3, :));
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(k_in.wt));
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(k_in.bfa1bub2));
violins = violinplot(xConc, gConc, 'Bandwidth', 1, 'Width', 0.4, 'ShowData', false);
violins(2).ViolinColor = colors(2, :);
violins(3).ViolinColor = colors(1, :);
violins(4).ViolinColor = colors(2, :);

plot([0.6 2.4], log10([1 1] .* pMin(2)), '--k');
plot([0.6 2.4], log10([1 1] .* pMax(2)), '--k');

plot([2.6 4.4], log10([1 1] .* pMin(3)), '--k');
plot([2.6 4.4], log10([1 1] .* pMax(3)), '--k');

xlim([0.5 4.5]);
ylim([-5.5 2.5]);
ylabel('nM Kip2^{-1} s^{-1}');

currAx = gca();
currAx.XTickLabel = {'WT', '\it{bfa1{\Delta}bub2{\Delta}}', 'WT', '\it{bfa1{\Delta}bub2{\Delta}}'};
currAx.YTick = -5:5;
currAx.YTickLabel = {'10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{-0}', '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}'};
currAx.XLabel.Color = [0 0 0];
currAx.YLabel.Color = [0 0 0];
currAx.XColor = [0 0 0];
currAx.YColor = [0 0 0];
xtickangle(45)

legend([violins(1).ViolinPlot, violins(2).ViolinPlot], {'WT', '\it{bfa1{\Delta}bub2{\Delta}}'}, 'Location', 'NorthWest');
legend('boxoff');

applyPaperFormatting

currFig = gcf();
currFig.Position(3) = 0.7*currFig.Position(3);
currFig.Color = 'white';

fprintf('Mean k_in: %g,\tStd: %g\n', mean(k_in.wt), std(k_in.wt));
fprintf('Median k_in: %g\n', median(k_in.wt));

fprintf('Mean k_in bfa1bub2del: %g,\tStd: %g\n', mean(k_in.bfa1bub2), std(k_in.bfa1bub2));
fprintf('Median k_in bfa1bub2del: %g\n', median(k_in.bfa1bub2));

if exportPlots
    print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_kon_kin.pdf'], '-dpdf');
    export_fig([parameterEstimateFigureResultFolder filesep 'png' filesep prefix '_violin_kon_kin.png'], '-q101', '-r300');
end

%%

r_on.wt   = (wtData.plotSamples(2, :) .* wtData.plotSamples(1, :) ./ 140 .* 61);
r_on.bfa1bub2 = (bfa1bub2Data.plotSamples(2, :) .* bfa1bub2Data.plotSamples(1, :) ./ 140 .* 61);
figure;

xConc = [];
gConc = [];
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(r_on.wt));
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(r_on.bfa1bub2));

fprintf('Mean r_on: %g,\t%g\n', mean(r_on.wt), std(r_on.wt));
fprintf('Median r_on: %g\n', median(r_on.wt));

fprintf('Mean r_on bfa1bub2: %g\t%g\n', mean(r_on.bfa1bub2), std(r_on.bfa1bub2));
fprintf('Median r_on bfa1bub2: %g\n', median(r_on.bfa1bub2));

r_in.wt   = (wtData.plotSamples(3, :) .* wtData.plotSamples(1, :) ./ 140 .* 61);
r_in.bfa1bub2 = (bfa1bub2Data.plotSamples(3, :) .* bfa1bub2Data.plotSamples(1, :) ./ 140 .* 61);

[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(r_in.wt));
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(r_in.bfa1bub2));
violins = violinplot(xConc, gConc, 'Bandwidth', 1, 'Width', 0.4, 'ShowData', false);
violins(1).ViolinColor = colors(1, :);
violins(2).ViolinColor = colors(2, :);
violins(3).ViolinColor = colors(1, :);
violins(4).ViolinColor = colors(2, :);

plot([0.6 2.4], log10([1 1] .* pMin(2) .* pMin(1) ./ 140 .* 61), '--k');
plot([0.6 2.4], log10([1 1] .* pMax(2) .* pMax(1) ./ 140 .* 61), '--k');

plot([2.6 4.4], log10([1 1] .* pMin(3) .* pMin(1) ./ 140 .* 61), '--k');
plot([2.6 4.4], log10([1 1] .* pMax(3) .* pMax(1) ./ 140 .* 61), '--k');

xlim([0.5 4.5]);
ylim([-3.5 4.5]);
ylabel('s^{-1}');

legend([violins(1).ViolinPlot, violins(2).ViolinPlot], {'WT', 'Kip2-bfa1bub2'}, 'Location', 'NorthWest');
legend('boxoff');

currAx = gca();
currAx.XTickLabel = {'WT', 'Kip2-bfa1bub2'};
currAx.YTick = -4:6;
currAx.YTickLabel = {'10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{-0}', '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}'};
currAx.XLabel.Color = [0 0 0];
currAx.YLabel.Color = [0 0 0];
currAx.XColor = [0 0 0];
currAx.YColor = [0 0 0];
xtickangle(45)
applyPaperFormatting

currFig = gcf();
currFig.Position(3) = 0.7*currFig.Position(3);
currFig.Color = 'white';

fprintf('Mean r_in: %g,\tStd: %g\n', mean(k_in.wt), std(k_in.wt));
fprintf('Median r_in: %g\n', median(k_in.wt));

fprintf('Mean r_in bfa1bub2: %g,\tStd: %g\n', mean(k_in.bfa1bub2), std(k_in.bfa1bub2));
fprintf('Median r_in bfa1bub2: %g\n', median(k_in.bfa1bub2));

if exportPlots
    print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_ron_rin.pdf'], '-dpdf');
    export_fig([parameterEstimateFigureResultFolder filesep 'png' filesep prefix '_violin_ron_rin.png'], '-q101', '-r300');
end
%%
k_out.wt   = (wtData.plotSamples(4, :));
k_out.bfa1bub2 = (bfa1bub2Data.plotSamples(4, :));
figure;
xConc = [];
gConc = [];
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, k_out.wt);
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, k_out.bfa1bub2);

violins = violinplot(xConc, gConc, 'Bandwidth', 0.5, 'Width', 0.4, 'ShowData', false);
violins(2).ViolinColor = colors(2, :);

plot([0.6 2.4], [1 1] .* pMin(4), '--k');
plot([0.6 2.4], [1 1] .* pMax(4), '--k');

xlim([0.5 2.5]);
ylim([0.5 14]);
ylabel('k_{out} (s^{-1})');


currAx = gca();
currAx.XTickLabel = {'WT', '\it{bfa1{\Delta}bub2{\Delta}}'};
currAx.XLabel.Color = [0 0 0];
currAx.YLabel.Color = [0 0 0];
currAx.XColor = [0 0 0];
currAx.YColor = [0 0 0];
xtickangle(45)

applyPaperFormatting

currFig = gcf();
currFig.Position(3) = 0.4*currFig.Position(3);
currFig.Color = 'white';

fprintf('Mean k_out: %g\n', mean(k_out.wt));
fprintf('Median k_out: %g\n', median(k_out.wt));    

fprintf('Mean k_out bfa1bub2del: %g\n', mean(k_out.bfa1bub2));
fprintf('Median k_out bfa1bub2del: %g\n', median(k_out.bfa1bub2));    

if exportPlots
    print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_kout.pdf'], '-dpdf');
    export_fig([parameterEstimateFigureResultFolder filesep 'png' filesep prefix '_violin_kout.png'], '-q101', '-r300');
end