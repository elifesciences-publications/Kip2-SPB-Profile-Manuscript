function useSimulationForPrediction()
exportFigures = true;
%% Configure folders & data file locations
simResultFolder = ['..' filesep 'simulationResults'];
simResultFolderSep = [simResultFolder filesep];

fileName = 'simResults-20180706T153609.mat'; % Produced using Fig2B mode in "runSampling.m"
analyzedFileName = ['analyzed-' fileName];

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
optimizationWorkspace = load([simResultFolderSep fileName]);
simResults = load([simResultFolderSep analyzedFileName]);
nSims = size(simResults.results, 2);
% Re-enable warnings
warning('on'); 

%% Load experimental data
expResults = load([expResultFolderSep expResultFile]);

%% List experimental conditions
displayConditions(expResults.conditionResults);

%% Configure experimental conditions to compare against
expResultsToCompare = 9:14; % Experimental conditions to plot
expResultFit = 11; % Condition that was fit originally
dataLineWidth = 1.181;

%% Plot models against experimental data
modelB = 10200;
ylims = [-1e3 8e3];
figure;
hold on;
colorOrder = get(gca, 'ColorOrder');
for expIndex = expResultsToCompare
    c = expResults.conditionResults{expIndex};
    fprintf('\nPlotting experimental condition: %s', c.condition);

    % SPB aligned
    shiftedLocations        = c.lengthVectorCenterShifted;
    SPBlocation             = (c.meanMeanOffset-0.5)*(4/30);
    shiftedLocations        = -(shiftedLocations - SPBlocation);
    greenIntensitiesShifted = c.greenIntensitiesCenterShifted;
    
    modelA = mean(c.backgroundIntensities);
    
    greenData = greenIntensitiesShifted;

    meanGreenData = mean(greenData, 2, 'omitnan');
    stdGreenData  = std(greenData, 0, 2, 'omitnan');
    nMTs = sum(~isnan(greenData),2);
    semGreenData = stdGreenData./sqrt(nMTs);
    
    [zero, zeroLoc] = min(abs(shiftedLocations));
    if zero > 1e-3
        error('SPB location failure');
    end

    shadedErrorBar(shiftedLocations, meanGreenData - modelA, semGreenData, 'lineProps', {'Color', 0.5*[1 1 1], 'LineWidth', dataLineWidth});                

    fitDataIndices = (shiftedLocations >= -0.07 & shiftedLocations < (shiftedLocations(zeroLoc - 2*c.meanMeanOffset) + 0.07));
    fitDataXcoords = shiftedLocations(fitDataIndices);
    
    dataPeakToPeakLengthInSites = round(round(mean(unique(c.lengthsUM))/0.133)*0.133/0.008);
    
    found = false;
    for simIndex = 1:nSims
        result = simResults.results{simIndex};
        p = optimizationWorkspace.parameterCombinations(:, simIndex);
        if p(1) == dataPeakToPeakLengthInSites
            xCoords = result.xCoords; 
            simMean = result.weightedIntensities;
            modelMeanFluo = modelB*simMean;            
            modelAtData = interp1(xCoords, modelMeanFluo, fitDataXcoords);
            
            if any(isnan(modelAtData))
                error('Model / data mismatch!');
            end
            
            if expIndex == expResultFit
                color = colorOrder(2, :);
            else
                color = colorOrder(1, :);
            end
            
            plot(fitDataXcoords, modelAtData, 'LineWidth', dataLineWidth, 'Color', color);
            
            leftOuterIndices =  xCoords <= min(fitDataXcoords);
            rightOuterIndices = xCoords >= max(fitDataXcoords);
            
            plot(xCoords(leftOuterIndices), modelMeanFluo(leftOuterIndices),  '--', 'LineWidth', dataLineWidth, 'Color', color);
            plot(xCoords(rightOuterIndices), modelMeanFluo(rightOuterIndices), '--', 'LineWidth', dataLineWidth, 'Color', color);
            found = true;
            break;
        end
    end
    if ~found
        error('No corresponding model found');
    end
    
end
fprintf('\n');
xlims = [-0.5, 2.8];
xlim(xlims);
ylim(ylims);
plot([0 0], ylims, '--', 'Color', [1 1 1] * 0.5);
plot(xlims, [0 0], '--', 'Color', [1 1 1] * 0.5);  
plot(xlims, modelB*[1 1], '--', 'Color', [1 1 1] * 0.5);  
xlabel('Distance from SPB ({\mu}m)');
ylabel(sprintf('Kip2-3xsfGFP fluorescence\n along aMT (a.u.)'));
applyPaperFormatting();
currentFig = gcf();
currentFig.Position(3) = 1.25 * currentFig.Position(3);
currentFig.Position(4) = 0.8 * currentFig.Position(4);

currentAxis = gca;
currentAxis.YTick = 0:2000:8000;
currentAxis.XLabel.Color = [0 0 0];
currentAxis.YLabel.Color = [0 0 0];
currentAxis.XColor = [0 0 0];
currentAxis.YColor = [0 0 0];

if exportFigures
    print([figureDirectorySep 'pdf' filesep 'PredictedBins.pdf'], '-dpdf');
    export_fig([figureDirectorySep 'png' filesep 'PredictedBins.png'], '-q101', '-r300');
end