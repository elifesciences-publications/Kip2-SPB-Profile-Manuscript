clear
close all
clc

%%
myFolder = fileparts(mfilename('fullpath'));
%%
simFileNames = {
     'simResults-20190221T061021.mat'
};
simResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'simulationResults'];

debug = false;
exportFigures = false;
prefix = 'all';

%%
expResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'analyzedData'];
expResultFile = [expResultFolder filesep 'binnedProfiles.mat'];

%%
parameterEstimateResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'simulationResults' filesep 'parameterEstimates'];
parameterEstimateFigureResultFolder = [myFolder filesep '..' filesep '..' filesep '..' filesep 'figures' filesep 'model'];

%% Load experimental results
expResults = load(expResultFile);

displayConditions(expResults.conditionResults);
% wt:             45:50 64:69 328:333
% wt mom:         55:60       338:343
% S63A:                 80:85
% S63A mom:             90:95
% bfa1del:                    346:351
% bfa1del mom:                355:360
% bub2del:                    365:370
% bub2del mom:                375:380
% bfa1/bub2del:               385:390
% bfa1/bub2del mom:           394:399

% [64:69 80:85]

% 346:351 365:370 
% 80:85 385:3901
% 80:85 385:390
%64:69
%3*ones(1,6) 
% 
% 45:50 64:69
% 
expResultsToCompare = [ 45:50 80:85 385:390];%[ ];%[ 45:50]% 64:69 328:333 ];% 64:69 328:333] %[45:50];%[64:69];%[45:50 64:69]; %45:50
expResultsNormalizeFrom = [ 1*ones(1,6) 3*ones(1,6) 35*ones(1,6)];% 3*ones(1,6) 35*ones(1,6)];% 3*ones(1,6) 35*ones(1,6) ];
%
expResultsNormalizeToReference = 1*ones(1,length(expResultsNormalizeFrom));%zeros(size(expResultsToCompare)); %0.*[1 1 1 1 1 1];
expResultsPlotID = [];
for i = 1:(length(expResultsToCompare)/6)
    expResultsPlotID = [expResultsPlotID ([1 1 1 1 1 1]*(i*10))];
end

nExpResultsToCompare = length(expResultsToCompare);
modelBs = [ 7000:500:13000 13100:100:14000 14500:500:17000 ];
nModelB = length(modelBs);

ylims = [-1.5 12];
scaleFluo = 1e-3;


%%
analyzedFileNames = simFileNames;
for i = 1:length(simFileNames)
    analyzedFileNames{i} = ['analyzed-' analyzedFileNames{i}];
end

optimizationWorkspaces = {};
simResultsAll = {};
binsAll = [];

%% Load simulation results
for i = 1:length(simFileNames)
    optimizationWorkspaces{end+1} = load([simResultFolder filesep simFileNames{i}]);
    
    if optimizationWorkspaces{end}.runSingleDebug || optimizationWorkspaces{end}.runSingleParallel
        warning(['File "' simFileNames{i} '" contains incomplete parameter set since debug mode was active.']);
        if i == 1
            parameterCombinations = optimizationWorkspaces{i}.parameterCombinations(:, optimizationWorkspaces{end}.debugIndices);
        else
            parameterCombinations = [parameterCombinations optimizationWorkspaces{i}.parameterCombinations(:, optimizationWorkspaces{end}.debugIndices)];
        end
        simResults = load([simResultFolder filesep analyzedFileNames{i}]);
        simResultsAll = [simResultsAll simResults.results{optimizationWorkspaces{end}.debugIndices}];
    else
        if i == 1
            parameterCombinations = optimizationWorkspaces{i}.parameterCombinations;
        else
            parameterCombinations = [parameterCombinations optimizationWorkspaces{i}.parameterCombinations];
        end
        simResults = load([simResultFolder filesep analyzedFileNames{i}]);
        simResultsAll = [simResultsAll simResults.results];
    end
    
    binsAll = unique([binsAll; simResults.binsToCompute], 'rows');
end


%%
nParameterSets = size(parameterCombinations,2);
nBins = size(binsAll, 1);

%% Match simulated bins to experimental conditions
SSE_arr = NaN(nModelB, nExpResultsToCompare, nParameterSets );
dof_arr = NaN(nExpResultsToCompare, 1);

binIndexCache = NaN(1, nExpResultsToCompare);
if debug
    debugFig = figure();
    hold on;
    colors = get(gca,'colororder');
    debugFig2 = figure();
end

binIndex = 0;
for bin = binsAll'
    binIndex = binIndex + 1;
    for expIndex = 1:length(expResultsToCompare)
        c = expResults.conditionResults{expResultsToCompare(expIndex)};
        condition = c.condition;

        if isnan(binIndexCache(expIndex))
            [tokens, ~] = regexp(condition,'^(.+?)-bin-(.+?)-(.+?)$','tokens','match');
            lowerBinBound = str2double(tokens{1}{2});
            upperBinBound = str2double(tokens{1}{3});

            if sum(abs(bin - [lowerBinBound;upperBinBound])) < 0.05
                binIndexCache(expIndex) = binIndex;
            end
        end
    end
end

assert(all(isfinite(binIndexCache)));
%% Laser power adjustment
laserFitResults = {};
laserGofs = {};
laserFitParams = nan(2, nExpResultsToCompare);
for expIndex = 1:nExpResultsToCompare
    referenceCondition = expResultsNormalizeToReference(expIndex);
    fromCondition = expResultsNormalizeFrom(expIndex);
    
    if (referenceCondition > 0 && fromCondition > 0)
        cFrom = expResults.conditionResults{fromCondition};
        cRef = expResults.conditionResults{referenceCondition};
        [laserFitResult, laserGof] = laserPowerModel(cFrom, cRef);
        laserFitResults{expIndex} = laserFitResult;
        laserGofs{expIndex} = laserGof;
        laserFitParams(:, expIndex) = [laserFitResult.p1; laserFitResult.p2];
    end
end
laserFitParams

%%$

%%
tic
parfor expIndex = 1:nExpResultsToCompare
    for currentParameterSetIndex = 1:nParameterSets
        currentResult = simResultsAll{currentParameterSetIndex};
        currentParameters = parameterCombinations(:, currentParameterSetIndex);
        binIndex = binIndexCache(expIndex);
        
        currentBinCoords = currentResult.vmBin{binIndex}.meanAllAlignedBinnedProfileCoords;
        currentBinProfile = currentResult.vmBin{binIndex}.meanAllAlignedBinnedProfile;
        
        if all(isnan(currentBinCoords))
            continue;
        end
        
        if expResultsNormalizeToReference(expIndex) > 0 && expResultsNormalizeFrom(expIndex) > 0
            laserFitParamsLocal = laserFitParams(:, expIndex);
        else
            laserFitParamsLocal = [];
        end
        
        c = expResults.conditionResults{expResultsToCompare(expIndex)};
        
        [fitDataXcoords, modelAtData, fitDataMean, fitDataSEM] = computeModelAtData(currentBinCoords, currentBinProfile, laserFitParamsLocal, c );

        for modelBindex = 1:nModelB             
            modelAtDataFluo = modelBs(modelBindex) * modelAtData;
            if debug
                figure(debugFig2);
                clf;
                hold on;
                plot(shiftedLocations, meanGreenData);
                plot(fitDataXcoords, fitDataMean, 'r');
                plot(fitDataXcoords, modelAtDataFluo, 'k');
            end
            SSE_sum = sum((modelAtDataFluo - fitDataMean).^2 ./ fitDataSEM.^2);
            
            if isnan(dof_arr(expIndex))
                dof_arr(expIndex) = length(modelAtDataFluo)
            elseif dof_arr(expIndex) ~= length(modelAtDataFluo)
                error('inconsistent DoFs');
            end

            SSE_arr(modelBindex, expIndex, currentParameterSetIndex) = SSE_sum;
        end
    end
end

toc

%%
totalDofs = sum(dof_arr);
[availableModelBs, maxSSE, modelBsse, bestSSE, I] = plotProfilePlots(prefix, SSE_arr, totalDofs, nExpResultsToCompare, parameterCombinations, modelBs, exportFigures, parameterEstimateFigureResultFolder);
%%
[~, bestModelBoverallIndex] = min(modelBsse);
[bestSSE2, bestModelBIndex] = min(bestSSE, [], 1);
[bestSSE3, bestParameterIndices] = min(bestSSE, [], 2);

viableParameterSetThresholds = chi2inv(0.95, dof_arr - 5)';
viableModelBs = bestSSE < viableParameterSetThresholds;
%%
parameterEstimatesIndividual = table();
compatibleParameterSetsInExpIndividual = cell(nExpResultsToCompare, 1);

for expIndex = 1:length(expResultsToCompare)
    currentSSEs = squeeze(SSE_arr(:, expIndex, :));
    currentMaxSSE = viableParameterSetThresholds(expIndex);
    
    [compatibleModelBs, compatibleParameterSets] = find(currentSSEs <= currentMaxSSE); 
    uniqueCompatibleModelBs = unique(compatibleModelBs);
    uniqueCompatibleParameterSets = unique(compatibleParameterSets);
    
    c = expResults.conditionResults{expResultsToCompare(expIndex)};
    condition = c.condition;
    
    bestParameterSetIndex = I(bestModelBIndex(expIndex), expIndex);
    binIndex = binIndexCache(expIndex);
    
    currentResult     = simResultsAll{bestParameterSetIndex};
    currentParameters = parameterCombinations(:, bestParameterSetIndex);

    fprintf('\n\nExperiment %i (condition %3i)\n', expIndex, expResultsToCompare(expIndex));
    fprintf('----------------------------\n');
    fprintf('Best modelB: %g (index %i)\n', modelBs(bestModelBIndex(expIndex)), bestModelBIndex(expIndex));
    fprintf('Best parameter set index: %i\n', bestParameterSetIndex);
    fprintf('Kip2 molecules: %i \n', currentParameters(1));
    fprintf('k_on: %g (s^-1 nM^-1)\n', currentParameters(2));
    fprintf('k_in: %g (s^-1 nM^-1)\n', currentParameters(3));
    fprintf('k_out: %g (s^-1)\n', currentParameters(4));
    
    currentBinCoords  = currentResult.vmBin{binIndex}.meanAllAlignedBinnedProfileCoords;
    currentBinProfile = currentResult.vmBin{binIndex}.meanAllAlignedBinnedProfile;
    
    if expResultsNormalizeToReference(expIndex) > 0 && expResultsNormalizeFrom(expIndex) > 0
        laserFitParamsLocal = laserFitParams(:, expIndex);
    else
        laserFitParamsLocal = [];
    end   
    [fitDataXcoords, modelAtData, fitDataMean, fitDataSEM, currentBinCoords, currentBinProfile, shiftedLocations,  meanGreenData, semGreenData] = computeModelAtData(currentBinCoords, currentBinProfile, laserFitParamsLocal, c );
    
    modelAtDataFluo = modelBs(bestModelBIndex(expIndex)) * modelAtData;
    modelAtDataFluoExtended = modelBs(bestModelBIndex(expIndex)) * currentBinProfile;

    figure(expResultsPlotID(expIndex));
    
    subplot(7,1,[1 2]);
    hold on;
    colors = get(gca,'colororder');
    
    H = shadedErrorBar(shiftedLocations, meanGreenData .* scaleFluo, semGreenData .* scaleFluo, 'lineProps', {'k', 'LineWidth', 1.5'});  
    %hAllData = plot(shiftedLocations, meanGreenData, 'Color', [1 1 1] .* 0.8, 'LineWidth', 1.5);
    %hFitData = plot(fitDataXcoords, fitDataMean, 'Color', [1 1 1] .* 0.5, 'LineWidth', 1.5);
    hModel = plot(fitDataXcoords, modelAtDataFluo .* scaleFluo, 'r', 'LineWidth', 1.5);
    
    
    modelIndices = currentBinCoords < fitDataXcoords(1) & currentBinCoords > shiftedLocations(end);
    hModel2 = plot(currentBinCoords(modelIndices), modelAtDataFluoExtended(modelIndices) .* scaleFluo, 'r--', 'LineWidth', 1.5);
    modelIndices = currentBinCoords > fitDataXcoords(end) & currentBinCoords < shiftedLocations(1);
    hModel3 = plot(currentBinCoords(modelIndices), modelAtDataFluoExtended(modelIndices) .* scaleFluo, 'r--', 'LineWidth', 1.5);
    
    p = currentParameters;
    legend([H.mainLine, hModel], {condition, 'Model'}, 'Location', 'NorthWest', 'AutoUpdate', 'off');

    %legend([H.mainLine, hModel], {condition, sprintf('nKip2Free = %g molecules\nk_{on} = %g\nk_{in} = %g\nk_{out} = %g\nB = %g', p(1:4),  modelBs(bestModelBIndex(expIndex)))}, 'Location', 'NorthWest');
    legend boxoff
    ylim(ylims);
    xlim([-0.5 4]);
    xlabel('Distance from SPB ({\mu}m)');
    ylabel('Occupancy');    
    
    x = mean(c.lengthsUM);
    
    subplot(7,1,3); % B
    hold on;
    plot(x, modelBs(bestModelBIndex(expIndex)), 'xk', 'LineWidth', 0.8);
    
    if ~isempty(uniqueCompatibleModelBs)
        y = modelBs(bestModelBIndex(expIndex)); 
        neg = y - modelBs(uniqueCompatibleModelBs(1));
        pos = modelBs(uniqueCompatibleModelBs(end)) - y;
        errorbar(x,y,neg,pos,'k','LineWidth', 0.8);
    end
    
    ylabel('B (a.u.)');
    xlim([-0.5 4]);
    
    subplot(7,1,4);
    hold on;
    compatibleParameterSetsInExpIndividual{expIndex} = compatibleParameterSets;

    plotEstimateErrorBar(parameterCombinations, uniqueCompatibleParameterSets, p, 1, x, 'k');
    ylabel('Kip2 (# molecules)');
    xlim([-0.5 4]);
    
    subplot(7,1,5);
    hold on;
    y = p(2)*p(1); 
    plot(x, y, 'xr', 'LineWidth', 0.8);
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(2, uniqueCompatibleParameterSets) .* parameterCombinations(1, uniqueCompatibleParameterSets));
        
        neg = y - uniqueParameterValues(1);
        pos = uniqueParameterValues(end) - y;
        errorbar(x,y,neg,pos,'r','LineWidth', 0.8);
    end
    
    y = p(3)*p(1); 
    plot(x, y, 'xb', 'LineWidth', 0.8);
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(3, uniqueCompatibleParameterSets) .* parameterCombinations(1, uniqueCompatibleParameterSets));
        
        neg = y - uniqueParameterValues(1);
        pos = uniqueParameterValues(end) - y;
        errorbar(x,y,neg,pos,'b','LineWidth', 0.8);
    end
    
    %plotEstimateErrorBar(parameterCombinations, uniqueCompatibleParameterSets, p, 2, x);
    ylabel('r_{max} (s^{-1})');
    xlim([-0.5 4]);
    
    subplot(7,1,6);
    hold on;
    plotEstimateErrorBar(parameterCombinations, uniqueCompatibleParameterSets, p, 2, x, 'r');
    plotEstimateErrorBar(parameterCombinations, uniqueCompatibleParameterSets, p, 3, x, 'b');
    ylabel('k_{in}, k_{on} (TODO)');
    xlim([-0.5 4]);
    
    % k_on, k_in
    subplot(7,1,7);
    hold on;
    plotEstimateErrorBar(parameterCombinations, uniqueCompatibleParameterSets, p, 4, x, 'k');    
    ylabel('k_{out} (s^{-1})');
    xlim([-0.5 4]);
    % k_out
    xlabel('Microtubule length ({\mu}m))');
    
    
    parameterEstimateIndividual = struct();
    parameterEstimateIndividual.conditionIndex = expResultsToCompare(expIndex);
    parameterEstimateIndividual.conditionName = {condition};
    parameterEstimateIndividual.meanLength = x;
    parameterEstimateIndividual.bestB = modelBs(bestModelBIndex(expIndex));
    parameterEstimateIndividual.B_CI = [NaN NaN];
    if ~isempty(uniqueCompatibleModelBs)
        parameterEstimateIndividual.B_CI = [modelBs(uniqueCompatibleModelBs(1)) modelBs(uniqueCompatibleModelBs(end))];
    end
    parameterEstimateIndividual.bestNKip2 = p(1); 
    parameterEstimateIndividual.nKip2CI = [NaN NaN];
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(1, uniqueCompatibleParameterSets));
        parameterEstimateIndividual.nKip2CI = [uniqueParameterValues(1) uniqueParameterValues(end)];
    end
    parameterEstimateIndividual.bestKOn = p(2); 
    parameterEstimateIndividual.kOnCI = [NaN NaN];
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(2, uniqueCompatibleParameterSets));
        parameterEstimateIndividual.kOnCI = [uniqueParameterValues(1) uniqueParameterValues(end)];
    end
    parameterEstimateIndividual.bestKIn = p(3); 
    parameterEstimateIndividual.kInCI = [NaN NaN];
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(3, uniqueCompatibleParameterSets));
        parameterEstimateIndividual.kInCI = [uniqueParameterValues(1) uniqueParameterValues(end)];
    end
    parameterEstimateIndividual.bestKOut = p(4); 
    parameterEstimateIndividual.kOutCI = [NaN NaN];
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(4, uniqueCompatibleParameterSets));
        parameterEstimateIndividual.kOutCI = [uniqueParameterValues(1) uniqueParameterValues(end)];
    end
    parameterEstimateIndividual.bestLoadingRatio = p(3)/p(2); 
    parameterEstimateIndividual.loadingRatioCI = [NaN NaN];
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(3, uniqueCompatibleParameterSets) ./ parameterCombinations(2, uniqueCompatibleParameterSets));
        parameterEstimateIndividual.loadingRatioCI = [uniqueParameterValues(1) uniqueParameterValues(end)];
    end
    parameterEstimateIndividual.bestROnMax = p(1)*p(2); 
    parameterEstimateIndividual.rOnMaxCI = [NaN NaN];
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(1, uniqueCompatibleParameterSets) .* parameterCombinations(2, uniqueCompatibleParameterSets));
        parameterEstimateIndividual.rOnMaxCI = [uniqueParameterValues(1) uniqueParameterValues(end)];
    end
    parameterEstimateIndividual.bestRInMax = p(1)*p(3); 
    parameterEstimateIndividual.rInMaxCI = [NaN NaN];
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(1, uniqueCompatibleParameterSets) .* parameterCombinations(3, uniqueCompatibleParameterSets));
        parameterEstimateIndividual.rInMaxCI = [uniqueParameterValues(1) uniqueParameterValues(end)];
    end
    
    
    parameterEstimatesIndividual = [parameterEstimatesIndividual; struct2table(parameterEstimateIndividual)];
end

%%
for i=unique(expResultsPlotID)
    
    figure(i);
    subplot(7,1,[1 2]);
    xlims = xlim();
    
    plot(xlims, [0 0], 'k--', 'HandleVisibility', 'off');
    plot([0 0], ylims, 'k--', 'HandleVisibility', 'off');
    
    subplot(7,1,3);
    plot(xlims, [1 1] .* min(modelBs), '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    plot(xlims, [1 1] .* max(modelBs), '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    ylim([min(modelBs)*0.5 2*max(modelBs)]);

    for j = 3:6
        subplot(7, 1, j);
        ax = gca();
        ax.YScale = 'log';
    end
    
    subplot(7, 1, 4);
    testedParameterRange = [parameterCombinations(1, :)];
    minTested = min(testedParameterRange);
    maxTested = max(testedParameterRange);
    plot(xlims, [1 1] .* minTested, '--k', 'LineWidth', 0.8);
    plot(xlims, [1 1] .* maxTested, '--k', 'LineWidth', 0.8);
    ylim([minTested*0.5 maxTested*2]);
    

    
    subplot(7, 1, 5);
    testedParameterRange = [(parameterCombinations(1, :) .* parameterCombinations(2, :)) (parameterCombinations(1, :) .* parameterCombinations(3, :))];
    minTested = min(testedParameterRange);
    maxTested = max(testedParameterRange);
    plot(xlims, [1 1] .* minTested, '--k', 'LineWidth', 0.8);
    plot(xlims, [1 1] .* maxTested, '--k', 'LineWidth', 0.8);
    ylim([minTested*0.5 maxTested*2]);
    
    subplot(7, 1, 6);
    testedParameterRange = [parameterCombinations(2, :)];
    minTested1 = min(testedParameterRange);
    maxTested1 = max(testedParameterRange);
    plot(xlims, [1 1] .* minTested1, '--r', 'LineWidth', 0.8);
    plot(xlims, [1 1] .* maxTested1, '--r', 'LineWidth', 0.8);
    
    testedParameterRange = [parameterCombinations(3, :)];
    minTested2 = min(testedParameterRange);
    maxTested2 = max(testedParameterRange);
    plot(xlims, [1 1] .* minTested2, '--b', 'LineWidth', 0.8);
    plot(xlims, [1 1] .* maxTested2, '--b', 'LineWidth', 0.8);
    
    ylim([min(minTested1,minTested2)*0.5 max(maxTested1,maxTested2)*2]);
    
    subplot(7,1,7);
    plot(xlims, [1 1] .* min(parameterCombinations(4, :)), '--k', 'LineWidth', 0.8);
    plot(xlims, [1 1] .* max(parameterCombinations(4, :)), '--k', 'LineWidth', 0.8);
    ylim([-1 max(parameterCombinations(4, :))*1.1]);
    
    applyPaperFormatting
    
    currFig = gcf();
    currFig.Position(4) = 2.5*currFig.Position(4);
    currFig.Color = 'white';
    movegui(currFig, 'center');
    drawnow();
    if exportFigures
        print([parameterEstimateFigureResultFolder filesep 'pdf' filesep 'paper-' prefix '-'  int2str(i) '.pdf'], '-dpdf');
        export_fig([parameterEstimateFigureResultFolder filesep 'png' filesep  'paper-' prefix '-'  int2str(i) '.png'], '-q101', '-r300');
	end
end

%% Plot wt beeswarm plots with sampling
pMin = min(parameterCombinations');
pMax = max(parameterCombinations');

nPlotSamples = 20000;
plotSamples = nan(4, nPlotSamples);
plotModelBs = nan(1, nPlotSamples);
sampleCondition = nan(1, nPlotSamples);

nProfilesPerCondition = nan(nExpResultsToCompare, 1);
for expIndex = 1:nExpResultsToCompare
    nProfilesPerCondition(expIndex) = size(expResults.conditionResults{expResultsToCompare(expIndex)}.lengths, 1);
end
profilesPerConditionTotal = nProfilesPerCondition ./ sum(nProfilesPerCondition);

likelihood_arr = exp(-0.5*SSE_arr);
likelihood_arr(isnan(likelihood_arr)) = 0;

parfor i = 1:nPlotSamples
    binFound = 0;
    sum_j = 0;
    prop = rand();
    for j = 1:nExpResultsToCompare
        jnext = j;
        sum_j1 = sum_j;
        sum_j = sum_j + profilesPerConditionTotal(j);
        if sum_j >= prop && prop > sum_j1
            binFound = j;
            break;
        end
    end
    if ~binFound
        error('No sample bin found?');
    end
    sampleCondition(i) = binFound;
    
    likelihood_arr_exp = squeeze(likelihood_arr(:, binFound, :));
    likelihood_arr_exp = likelihood_arr_exp ./ sum(likelihood_arr_exp(:));
    
    [modelBIndices, parameterIndices] = find(likelihood_arr_exp >= 0.1/nPlotSamples);
    likelihood_arr_exp = likelihood_arr_exp(likelihood_arr_exp >= 0.1/nPlotSamples);
    
    if sum(likelihood_arr_exp(:)) < 0.99
        error('Sampling error too large (> 1%)');
    end
    
    likelihood_arr_exp = likelihood_arr_exp ./ sum(likelihood_arr_exp);
    
    sampleFound = 0;
    nSamplesInLikelihoodArr = length(likelihood_arr_exp);
    
    sum_j = 0;
    prop = rand();
    for j = 1:nSamplesInLikelihoodArr
        jnext = j;
        sum_j1 = sum_j;
        sum_j = sum_j + likelihood_arr_exp(j);
        if sum_j >= prop && prop > sum_j1
            sampleFound = j;
            break;
        end
    end
    
    if ~sampleFound
        error('No sample bin found?');
    end
    
    plotSamples(:, i) = parameterCombinations(:, parameterIndices(sampleFound)); 
    plotModelBs(:, i) = modelBs(modelBIndices(sampleFound));
end

plotViolinPlots([prefix '_sampled'], plotSamples, pMin, pMax, exportFigures, parameterEstimateFigureResultFolder);

save([parameterEstimateResultFolder filesep prefix '_params.mat'], 'prefix', 'plotSamples', 'parameterCombinations');