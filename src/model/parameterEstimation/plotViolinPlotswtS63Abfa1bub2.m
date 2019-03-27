clear
close all
clc

myFolder = fileparts(mfilename('fullpath'));

exportPlots = true;

parameterEstimateResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'simulationResults' filesep 'parameterEstimates'];
parameterEstimateFigureResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'figures' filesep 'model'];

wtData = load([parameterEstimateResultFolder filesep 'wt_params.mat']);
S63AData = load([parameterEstimateResultFolder filesep 'S63A_params.mat']);
bfa1bub2Data = load([parameterEstimateResultFolder filesep 'bfa1bub2_params.mat']);

prefix = 'comparison-all';

pMin = min(wtData.parameterCombinations');
pMax = max(wtData.parameterCombinations');

%% Plot beeswarm plots
figure;
currAx = gca();
colors = currAx.ColorOrder;

concentrations.wt   = (wtData.plotSamples(1, :) ./ 140 .* 61);
concentrations.S63A = (S63AData.plotSamples(1, :) ./ 140 .* 61);
concentrations.bfa1bub2 = (bfa1bub2Data.plotSamples(1, :) ./ 140 .* 61);
[xConc, gConc] = addBoxPlotGroup([], [], concentrations.wt);
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, concentrations.S63A);
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, concentrations.bfa1bub2);

violins = violinplot(xConc, gConc, 'Bandwidth', 30, 'Width', 0.4, 'ShowData', false);
violins(2).ViolinColor = colors(2, :);
violins(3).ViolinColor = colors(3, :);
xlim([0.5 3.5]);
ylabel('[Kip2]_{total} (nM)');
plot([0.6 3.4], [1 1] .* pMin(1) ./ 140 .* 61, '--k');
plot([0.6 3.4], [1 1] .* pMax(1) ./ 140 .* 61, '--k');
ylim([0 230]);

currAx.XTickLabel = {'WT', 'Kip2-S63A', '\it{bfa1{\Delta}bub2{\Delta}}'};
currAx.XLabel.Color = [0 0 0];
currAx.YLabel.Color = [0 0 0];
currAx.XColor = [0 0 0];
currAx.YColor = [0 0 0];
xtickangle(45)
applyPaperFormatting

currFig = gcf();
currFig.Position(3) = 0.55*currFig.Position(3);
currFig.Color = 'white';


fprintf('\nMean [Kip2]: %g,\tStd:%g\n', mean(concentrations.wt), std(concentrations.wt));
fprintf('Median [Kip2]: %g\n', median(concentrations.wt));

fprintf('Mean [Kip2-S63A]: %g,\tStd:%g\n', mean(concentrations.S63A), std(concentrations.S63A));
fprintf('Median [Kip2-S63A]: %g\n', median(concentrations.S63A));

fprintf('Mean [Kip2 bfa1bub2]: %g,\tStd:%g\n', mean(concentrations.bfa1bub2), std(concentrations.bfa1bub2));
fprintf('Median [Kip2 bfa1bub2]: %g\n', median(concentrations.bfa1bub2));

if exportPlots
    print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_conc.pdf'], '-dpdf');
    export_fig([parameterEstimateFigureResultFolder filesep 'png' filesep prefix '_violin_conc.png'], '-q101', '-r300');
end

%%
k_out.wt   = (wtData.plotSamples(4, :));
k_out.S63A = (S63AData.plotSamples(4, :));
k_out.bfa1bub2 = (bfa1bub2Data.plotSamples(4, :));

figure;
xConc = [];
gConc = [];
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, k_out.wt);
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, k_out.S63A);
[xConc, gConc] = addBoxPlotGroup(xConc, gConc, k_out.bfa1bub2);

violins = violinplot(xConc, gConc, 'Bandwidth', 0.5, 'Width', 0.4, 'ShowData', false);
violins(2).ViolinColor = colors(2, :);
violins(3).ViolinColor = colors(3, :);


plot([0.6 3.4], [1 1] .* pMin(4), '--k');
plot([0.6 3.4], [1 1] .* pMax(4), '--k');

xlim([0.5 3.5]);
ylim([0.5 14]);
ylabel('k_{out} (s^{-1})');

xtickangle(45)
currAx = gca();
currAx.XTickLabel = {'WT', 'Kip2-S63A', '\it{bfa1{\Delta}bub2{\Delta}}'};
currAx.XLabel.Color = [0 0 0];
currAx.YLabel.Color = [0 0 0];
currAx.XColor = [0 0 0];
currAx.YColor = [0 0 0];


applyPaperFormatting

currFig = gcf();
currFig.Position(3) = 0.55*currFig.Position(3);
currFig.Color = 'white';

fprintf('\nMean k_out: %g\n', mean(k_out.wt));
fprintf('Median k_out: %g\n', median(k_out.wt));    

fprintf('Mean k_out S63A: %g\n', mean(k_out.S63A));
fprintf('Median k_out S63A: %g\n', median(k_out.S63A));    

fprintf('Mean k_out bfa1bub2: %g\n', mean(k_out.bfa1bub2));
fprintf('Median k_out bfa1bub2: %g\n', median(k_out.bfa1bub2)); 

if exportPlots
    print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_kout.pdf'], '-dpdf');
    export_fig([parameterEstimateFigureResultFolder filesep 'png' filesep prefix '_violin_kout.png'], '-q101', '-r300');
end
%%
k_on.wt       = (wtData.plotSamples(2, :));
k_on.S63A     = (S63AData.plotSamples(2, :));
k_on.bfa1bub2 = (bfa1bub2Data.plotSamples(2, :));

k_in.wt       = (wtData.plotSamples(3, :));
k_in.S63A     = (S63AData.plotSamples(3, :));
k_in.bfa1bub2 = (bfa1bub2Data.plotSamples(3, :));
%% Check for concentration differences
fprintf('\nConcentration difference: p = %g \t(H0: wt < S63A)\n', mean(concentrations.wt > concentrations.S63A));
fprintf('Concentration difference: p = %g \t(H0: wt < bfa1bub2)\n', mean(concentrations.wt > concentrations.bfa1bub2));

%% WT in vs on rate
fprintf('\n');
prob = mean(k_on.wt > k_in.wt);
if prob == 0
    fprintf('WT in rate > on rate: p < %g \t(H0: k_in < k_on)\n', 1/length(k_on.wt));
else
    fprintf('WT in rate > on rate: p = %g \t(H0: k_in < k_on)\n', prob);
end

%% S63A vs WT in / on rates
S63A_concentrations = unique(S63AData.plotSamples(1, :));
wt_concentrations = unique(wtData.plotSamples(1, :));
overlappingSamples = intersect(wt_concentrations, S63A_concentrations);
pValuesKon = zeros(size(overlappingSamples));
pValuesKin = zeros(size(overlappingSamples));
pValuesKout = zeros(size(overlappingSamples));
nSamples = zeros(size(overlappingSamples));

fprintf('\n');
i = 1;
for conc = overlappingSamples
    wt_samples   = find(wtData.plotSamples(1, :) == conc);
    S63A_samples = find(S63AData.plotSamples(1, :) == conc);
    
    wt_samples_rand   = wt_samples(randi(length(wt_samples), 1, 50000));
    S63A_samples_rand = S63A_samples(randi(length(S63A_samples), 1, 50000));
    
    pValuesKon(i) = mean(wtData.plotSamples(2, wt_samples_rand) >= S63AData.plotSamples(2, S63A_samples_rand));
    pValuesKin(i) = mean(wtData.plotSamples(3, wt_samples_rand) <= S63AData.plotSamples(3, S63A_samples_rand));
    nSamples(i) = min(length(wt_samples),length(S63A_samples));
    i = i + 1;
end
sharedConcentrations = overlappingSamples ./ 140 .* 61;
nSamples
pValuesKon
pValuesKin

if sum(pValuesKon(2:end)) == 0
    fprintf('S63A on rate > WT on rate (for concentrations >= %2.1f nM): p < %g \t(H0: S63A on rate <= WT on rate)\n',sharedConcentrations(2), 1./sum(nSamples(2:end)));
else
    fprintf('S63A on rate > WT on rate (for concentrations >= %2.1f nM): p = %g \t(H0: S63A on rate <= WT on rate)\n',sharedConcentrations(2), sum(pValuesKon(2:end) .* nSamples(2:end))/sum(nSamples(2:end)));
end

%%
bfa1bub2_concentrations = unique(bfa1bub2Data.plotSamples(1, :));
wt_concentrations = unique(wtData.plotSamples(1, :));
overlappingSamples = intersect(wt_concentrations, bfa1bub2_concentrations);
pValuesKon2 = zeros(size(overlappingSamples));
pValuesKin2 = zeros(size(overlappingSamples));
nSamples = zeros(size(overlappingSamples));

i = 1;
for conc = overlappingSamples
    wt_samples   = find(wtData.plotSamples(1, :) == conc);
    bfa1bub2_samples = find(bfa1bub2Data.plotSamples(1, :) == conc);
    
    wt_samples_rand   = wt_samples(randi(length(wt_samples), 1, 50000));
    bfa1bub2_samples_rand = bfa1bub2_samples(randi(length(bfa1bub2_samples), 1, 50000));
    
    pValuesKon2(i) = mean(wtData.plotSamples(2, wt_samples_rand) >= bfa1bub2Data.plotSamples(2, bfa1bub2_samples_rand));
    pValuesKin2(i) = mean(wtData.plotSamples(3, wt_samples_rand) <= bfa1bub2Data.plotSamples(3, bfa1bub2_samples_rand));
    nSamples(i) = min(length(wt_samples),length(bfa1bub2_samples));
    i = i + 1;
end
sharedConcentrations = overlappingSamples ./ 140 .* 61;
nSamples
pValuesKon2
pValuesKin2

if sum(pValuesKin2(2:end)) == 0
    fprintf('bfa1bub2del in rate < WT in rate (for concentrations >= %2.1f nM): p < %g \t(H0: bfa1bub2del in rate >= WT in rate)\n',sharedConcentrations(2), 1./sum(nSamples(2:end)));
else
    fprintf('bfa1bub2del in rate < WT in rate (for concentrations >= %2.1f nM): p = %g \t(H0: bfa1bub2del in rate >= WT in rate)\n',sharedConcentrations(2), sum(pValuesKin2(2:end) .* nSamples(2:end))/sum(nSamples(2:end)));
    
end
