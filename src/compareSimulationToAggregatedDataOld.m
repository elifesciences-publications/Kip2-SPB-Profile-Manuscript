function compareSimulationToAggregatedData()
% compareSimulationToAggregatedData: Uses binned profile data to
% compute the likelihood of model profiles, and creates the plots
% Fig2C and S8 in the manuscript.
%
%   Usage:
%   compareSimulationToAggregatedData()

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

exportFigures = true;
%% Configure folders & data file locations
simResultFolder = ['..' filesep 'simulationResults'];
simResultFolderSep = [simResultFolder filesep];

fileName1 =  'simResults-20180626T202346.mat'; % Produced using Fig2C-1 mode in "runSampling.m"
fileName2 = 'simResults-20180628T094544.mat';  % Produced using Fig2C-2 mode in "runSampling.m"

analyzedFileName1 = ['analyzed-' fileName1];
analyzedFileName2 = ['analyzed-' fileName2];

expResultFolder = ['..' filesep 'analyzedData'];
expResultFolderSep = [expResultFolder filesep];
expResultFile = 'binnedProfiles.mat';

figureDirectory = ['..' filesep 'figures' filesep 'model'];
figureDirectorySep = [figureDirectory filesep];
%% Create figure directories if necessary 
if exportFigures
    if ~exist(figureDirectory, 'dir')
        mkdir(figureDirectory);
    end

    if ~exist([figureDirectorySep 'png'], 'dir')
        mkdir([figureDirectorySep 'png']);
    end

    if ~exist([figureDirectorySep 'pdf'], 'dir')
        mkdir([figureDirectorySep 'pdf']);
    end
end
    

%% Load simulation data

% Since we are not on the cluster and the functions live in different locations, disable function handle warning
warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle'); 
optimizationWorkspace1 = load([simResultFolderSep fileName1]);
simResults1 = load([simResultFolderSep analyzedFileName1]);

optimizationWorkspace2 = load([simResultFolderSep fileName2]);
simResults2 = load([simResultFolderSep analyzedFileName2]);
% Re-enable warnings
warning('on'); 

parameterCombinations1 = optimizationWorkspace1.parameterCombinations;
parameterCombinations2 = optimizationWorkspace2.parameterCombinations;
parameterCombinations = [parameterCombinations1 parameterCombinations2];

% Convert molecules / PF to nM:
parameterCombinations(2, :) = parameterCombinations(2, :) ./ (140) .* 61;
simResults = struct;
simResults.results = [simResults1.results simResults2.results];

nSims = size(simResults.results, 2);

%% Load experimental data
expResults = load([expResultFolderSep expResultFile]);

%% List experimental conditions
displayConditions(expResults.conditionResults);

%% Configure experimental condition to compare against
expIndex = 11;

%% Compare models against experimental data
c = expResults.conditionResults{expIndex};
fprintf('\nSelected condition: %s\n', c.condition);

% Use SPB-aligned data
shiftedLocations        = c.lengthVectorCenterShifted;
SPBlocation             = (c.meanMeanOffset-0.5)*(4/30);
shiftedLocations        = -(shiftedLocations - SPBlocation);
greenIntensitiesShifted = c.greenIntensitiesCenterShifted;

modelA = mean(c.backgroundIntensities);
modelB = 10200;

greenData = greenIntensitiesShifted;
meanGreenData = mean(greenData, 2, 'omitnan');
stdGreenData  = std(greenData, 0, 2, 'omitnan');
nMTs = sum(~isnan(greenData),2);
semGreenData = stdGreenData./sqrt(nMTs);

xlimIndices = [find(any(isfinite(greenIntensitiesShifted),2),1) find(any(isfinite(greenIntensitiesShifted),2),1,'last')];
xlims = sort(shiftedLocations(xlimIndices));

[zero, zeroLoc] = min(abs(shiftedLocations));
if zero > 1e-3
    error('SPB location failure');
end


minSSE = Inf;
bestSim = 0;

SSE_arr = nan(1, nSims);
SSE_ndata = nan(1, nSims);
chi_arr = nan(1, nSims);


fitDataIndices = (shiftedLocations >= -0.07 & shiftedLocations < (shiftedLocations(zeroLoc - 2*c.meanMeanOffset) + 0.07));
fitDataXcoords = shiftedLocations(fitDataIndices);

fitDataMean = meanGreenData(fitDataIndices);
fitDataSEM  = semGreenData(fitDataIndices);   

for simIndex = 1:nSims
    result = simResults.results{simIndex};

    xCoords = result.xCoords; 
    simMean = result.weightedIntensities;
    modelMeanFluo = modelA + modelB*simMean;

    modelAtData = interp1(xCoords, modelMeanFluo, fitDataXcoords);

    if any(isnan(modelAtData))
        warning('Model data incompatiple with experimental profile');
        continue;
    end

    SSE_sum = sum((modelAtData - fitDataMean).^2 ./ fitDataSEM.^2);
    
    SSE_arr(simIndex) = SSE_sum;
    SSE_ndata(simIndex) = length(modelAtData);
    chi_arr(simIndex) = chi2inv(0.95, length(modelAtData) - 3); % + 1 for concentration comparison
    if ~isfinite(SSE_sum)
        error('Error: data mean contains NaN or data SEM is zero!');
    end

    if SSE_sum < minSSE
        minSSE = SSE_sum;
        bestSim = simIndex;
    end
end

if length(unique(chi_arr)) ~= 1
    error('Non-unique data length');
end

result = simResults.results{bestSim};
xCoords = result.xCoords; 
simMean = result.weightedIntensities;
modelMeanFluo = modelA + modelB*simMean;

modelAtData = interp1(xCoords, modelMeanFluo, fitDataXcoords);

%% Plot best model
baseline = modelA;
maxline = modelA + modelB;
ylims(1) = 0.9*baseline;
ylims(2) = 1.1*maxline;

figure();
hold on;

plot(shiftedLocations, meanGreenData);
shadedErrorBar(shiftedLocations, meanGreenData, semGreenData);                

plot([0 0], ylims, '--', 'Color', [1 1 1] * 0.5);
plot([1 1]*shiftedLocations(zeroLoc - 2*c.meanMeanOffset), ylims, '--', 'Color', [1 1 1] * 0.5);
plot(xlims, ([baseline baseline]), '--', 'Color', [1 1 1] * 0.5);
xlim(xlims);
ylim(ylims);
title(c.condition);
modelHandle = plot(fitDataXcoords, modelAtData, 'r');
leftOuterIndices =  xCoords < min(fitDataXcoords);
rightOuterIndices = xCoords > max(fitDataXcoords);
plot(xCoords(leftOuterIndices), modelMeanFluo(leftOuterIndices), 'r--');
plot(xCoords(rightOuterIndices), modelMeanFluo(rightOuterIndices), 'r--');
p = parameterCombinations(:, bestSim);
legend(modelHandle, sprintf('nKip2Free = %g nM\nk_{on} = %g nM^{-1} s^{-1} \nk_{in} = %g nM^{-1} s^{-1}\nk_{out} = %g s^{-1}\n', p(2:5)));
legend boxoff
applyPaperFormatting;
box off;
xlabel('Distance ({\mu}m)');
ylabel(sprintf('Kip2-3xsfGFP fluorescence\n along aMT (a.u.)'));
if exportFigures
    print([figureDirectorySep 'pdf' filesep 'FitToCentralBin.pdf'], '-dpdf');
    export_fig([figureDirectorySep 'png' filesep 'FitToCentralBin.png'], '-q101', '-r300');
end

%% Plot Likelihood Landscape
paramsOptimized = [2 4 5];
pMin = min(parameterCombinations,[],2);
pMax = max(parameterCombinations,[],2);
paramSpecs = readtable('paramSpecsConc.txt');
paramSpecs.bmin = pMin(paramsOptimized);
paramSpecs.bmax = pMax(paramsOptimized);
paramSpecs.bmin(1) = paramSpecs.bmin(1);
paramSpecs.bmax(1) = 600/140*61;
paramSpecs.bmax(2) = 2;
paramSpecs.bmin(3) = 2.5;
paramSpecs.bmax(3) = 5.5;

OutV = struct;
OutV.V = parameterCombinations(paramsOptimized, :)';
OutV.cost = SSE_arr';

plot_OutV2(OutV, paramSpecs, 'plotMode', 'sampleValues', 'nDoF', length(modelAtData), 'nParameters', 3, 'ticks', {50:100:250, [0.01 0.1 1], 3:5});
applyPaperFormatting
currFig = gcf();
currFig.Position(3:4) = round(0.9.*[1.67 1.5]  .* currFig.Position(3:4));
movegui(currFig, 'center');

if exportFigures
    export_fig([figureDirectorySep 'png' filesep 'LikelihoodLandscape.png'], '-q101', '-r300');
    currFig.PaperPositionMode = 'manual';
    orient(currFig, 'landscape')
    print([figureDirectorySep 'pdf' filesep 'LikelihoodLandscape.pdf'], '-dpdf');
end


%% Contrast maximum in rate vs out rate, and in rate constant vs total Kip2
pNew = OutV.V(:, 1) .* OutV.V(:, 2); % In rate in 1/s
figure;
cMap = parula();
cMap = flipud(cMap);

subplot(1,3,1);
colormap(cMap);
hold on;
scatter(pNew, OutV.V(:, 3), 30, 'filled', 'MarkerFaceAlpha', 0.05, 'MarkerFaceColor', [0.5 0.5 0.5]);

levelSets = [0.95, 0.75, 0.5, 0.25, 0.05];

minCost = min(OutV.cost);
maxCostValues = nan(size(levelSets));
for j = 1:length(levelSets)
    maxCostValues(j) = chi2inv(levelSets(j), length(modelAtData) - 3);
end
maxCost = max(maxCostValues);

for j = 1:length(levelSets)
    levelSetThreshold = maxCostValues(j);
    levelSetIndices = OutV.cost <= levelSetThreshold;
    if j < length(levelSets)
        levelSetIndices = levelSetIndices & (OutV.cost > maxCostValues(j + 1));
    end

    xSet = pNew(levelSetIndices);

    if ~isempty(xSet)
        ySet = OutV.V(levelSetIndices, 3);
        costSet = OutV.cost(levelSetIndices);
        normalizedCost = levelSetThreshold - minCost;
        normalizedCost = normalizedCost ./ (maxCost - minCost);

        scatter(xSet, ySet, 30, costSet, 'filled', 'MarkerFaceAlpha', (1 - normalizedCost) * 0.5 + 0.5);
    end
end
[~, maxLikelihoodIndex] = min(OutV.cost);
scatter(pNew(maxLikelihoodIndex), OutV.V(maxLikelihoodIndex, 3), 30, OutV.cost(maxLikelihoodIndex), 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 1);
plot([1 pNew(maxLikelihoodIndex)], OutV.V(maxLikelihoodIndex, 3) * [1 1], '--r');
plot(pNew(maxLikelihoodIndex) * [1 1], [1 OutV.V(maxLikelihoodIndex, 3)], '--r');

xlabel('k_{in} \cdot [Kip2]_{total} (s^{-1})');
ylabel('k_{out} (s^{-1})');

xlim([1 15]);
ylim([1 8]);
currentAxis = gca;
currentAxis.XLabel.Color = [0 0 0];
currentAxis.YLabel.Color = [0 0 0];
currentAxis.XColor = [0 0 0];
currentAxis.YColor = [0 0 0];    

subplot(1,3,2);
colormap(cMap);
hold on;
scatter(OutV.V(:, 1), 1./OutV.V(:, 2), 30, 'filled', 'MarkerFaceAlpha', 0.05, 'MarkerFaceColor', [0.5 0.5 0.5]);

for j = 1:length(levelSets)
    levelSetThreshold = maxCostValues(j);
    levelSetIndices = OutV.cost <= levelSetThreshold;
    if j < length(levelSets)
        levelSetIndices = levelSetIndices & (OutV.cost > maxCostValues(j + 1));
    end

    xSet = OutV.V(levelSetIndices, 1);

    if ~isempty(xSet)
        ySet = 1./OutV.V(levelSetIndices, 2);
        costSet = OutV.cost(levelSetIndices);
        normalizedCost = levelSetThreshold - minCost;
        normalizedCost = normalizedCost ./ (maxCost - minCost);

        scatter(xSet, ySet, 30, costSet, 'filled', 'MarkerFaceAlpha', (1 - normalizedCost) * 0.5 + 0.5);
    end
end

scatter(OutV.V(maxLikelihoodIndex, 1), 1/OutV.V(maxLikelihoodIndex, 2), 30, OutV.cost(maxLikelihoodIndex), 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 1);
plot([1 OutV.V(maxLikelihoodIndex, 1)], 1/OutV.V(maxLikelihoodIndex, 2) * [1 1], '--r');
plot(OutV.V(maxLikelihoodIndex, 1) * [1 1], [1 1/OutV.V(maxLikelihoodIndex, 2)], '--r');

currentAxis = gca;
currentAxis.XLabel.Color = [0 0 0];
currentAxis.YLabel.Color = [0 0 0];
currentAxis.XColor = [0 0 0];
currentAxis.YColor = [0 0 0];    

c = colorbar();
c.XTick = fliplr(maxCostValues);
c.Limits(2) = maxCost;
c.Color = [0 0 0];
for i = 1:length(maxCostValues)
    c.TickLabels{i} = sprintf('%g%%', 100*levelSets(length(maxCostValues) - i + 1));
end
c.Label.String = 'Confidence Region';

pos = get(c, 'Position');


set (c, 'Position', [pos(1)+0.1 pos(2) pos(3) pos(4)]);
xlabel('[Kip2]_{total} (nM)');
ylabel('k_{in}^{-1} (nM s)');
xlim([0 300]);
ylim([0 200]);


applyPaperFormatting;
currFig = gcf();
currFig.Position(3) = 2 .* currFig.Position(3);

movegui(currFig, 'center');
if exportFigures
    export_fig([figureDirectorySep 'png' filesep 'LikelihoodLandscapeMain.png'], '-q101', '-r300');
    currFig.PaperPositionMode = 'manual';
    orient(currFig, 'landscape')
    print([figureDirectorySep 'pdf' filesep 'LikelihoodLandscapeMain.pdf'], '-dpdf');
end