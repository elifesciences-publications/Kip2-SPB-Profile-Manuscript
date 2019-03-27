function [availableModelBs, maxSSE, modelBsse, bestSSE, I]  = plotProfilePlots(prefix, SSE_arr, totalDofs, nExpResultsToCompare, parameterCombinations, modelBs, exportFigures, parameterEstimateFigureResultFolder)
%% B consistency, wt
[bestSSE, I] = min(SSE_arr,[],3);

%SSE_arr2  = squeeze(sum(SSE_arr(, 2));

maxSSE = chi2inv(0.99, totalDofs - 1 - nExpResultsToCompare*4); % Estimate B over all data, and all other parameters per dataset
modelBsse = sum(bestSSE, 2);

figure();
semilogy(modelBs, modelBsse, 'k', 'LineWidth', 0.8);
hold on;
semilogy([min(modelBs) max(modelBs)], [1 1] .* maxSSE, 'r--', 'LineWidth', 0.8);
xlabel('B');
ylabel('Total SSE');
if length(modelBs) > 1
    xlim([min(modelBs) max(modelBs)]);
end

applyPaperFormatting
availableModelBs = modelBsse <= maxSSE;
availableModelBs = find(availableModelBs);

folderName = parameterEstimateFigureResultFolder;
figFileName = [prefix '_profile_modelB'];

if exportFigures
    print([folderName filesep 'pdf' filesep figFileName '.pdf'], '-dpdf');
    export_fig([folderName filesep 'png' filesep figFileName '.png'], '-q101', '-r300');
end
%% Kip2 concentration consistency
kip2Counts = unique(parameterCombinations(1, :));
nKip2Counts = length(kip2Counts);

nKip2SSEs = nan(nKip2Counts, 1);
for i = 1:nKip2Counts
    parameterIndices = parameterCombinations(1, :) == kip2Counts(i);
    
    minSSElocal = 0;
    for expIndex = 1:nExpResultsToCompare
        SSEs = squeeze(SSE_arr(:,expIndex,parameterIndices));
        minSSElocal = minSSElocal + min(min(SSEs));
    end
    nKip2SSEs(i) = minSSElocal;
end

figure();
semilogy(kip2Counts ./ 140 .* 61, nKip2SSEs, 'k', 'LineWidth', 0.8);
hold on;
semilogy([min(kip2Counts) max(kip2Counts)] ./ 140 .* 61, [1 1] .* maxSSE, 'r--', 'LineWidth', 0.8);
xlabel('[Kip2]_{total}');
ylabel('Total SSE');
xlim([min(kip2Counts) max(kip2Counts)] ./ 140 .* 61);
applyPaperFormatting

figFileName = [prefix '_profile_conc'];
if exportFigures
    print([folderName filesep 'pdf' filesep figFileName '.pdf'], '-dpdf');
    export_fig([folderName filesep 'png' filesep figFileName '.png'], '-q101', '-r300');
end

%% K on wt
konValues = unique(parameterCombinations(2, :));
nKonValues = length(konValues);

nKonSSEs = nan(nKonValues, 1);
for i = 1:nKonValues
    parameterIndices = parameterCombinations(2, :) == konValues(i);
    
    minSSElocal = 0;
    for expIndex = 1:nExpResultsToCompare
        SSEs = squeeze(SSE_arr(:,expIndex, parameterIndices));
        minSSElocal = minSSElocal + min(min(SSEs));
    end
    nKonSSEs(i) = minSSElocal;
end

figure();
loglog(konValues, nKonSSEs, 'k', 'LineWidth', 0.8);
hold on;
loglog([min(konValues) max(konValues)], [1 1] .* maxSSE, 'r--', 'LineWidth', 0.8);
xlabel('k_{on}');
ylabel('Total SSE');
xlim([min(konValues) max(konValues)]);
applyPaperFormatting

figFileName = [prefix '_profile_kon'];
if exportFigures
    print([folderName filesep 'pdf' filesep figFileName '.pdf'], '-dpdf');
    export_fig([folderName filesep 'png' filesep figFileName '.png'], '-q101', '-r300');
end


%% K in wt
kinValues = unique(parameterCombinations(3, :));
nkinValues = length(kinValues);

nKinSSEs = nan(nkinValues, 1);
for i = 1:nkinValues
    parameterIndices = parameterCombinations(3, :) == kinValues(i);
    
    minSSElocal = 0;
    for expIndex = 1:nExpResultsToCompare
        SSEs = squeeze(SSE_arr(:,expIndex,parameterIndices));
        minSSElocal = minSSElocal + min(min(SSEs));
    end
    nKinSSEs(i) = minSSElocal;
end
figure();
loglog(kinValues, nKinSSEs, 'k', 'LineWidth', 0.8);
hold on;
loglog([min(kinValues) max(kinValues)], [1 1] .* maxSSE, 'r--', 'LineWidth', 0.8);
xlabel('k_{in}');
ylabel('Total SSE');
xlim([min(kinValues) max(kinValues)]);
applyPaperFormatting

figFileName = [prefix '_profile_kin'];
if exportFigures
    print([folderName filesep 'pdf' filesep figFileName '.pdf'], '-dpdf');
    export_fig([folderName filesep 'png' filesep figFileName '.png'], '-q101', '-r300');
end

%% K out wt
koutValues = unique(parameterCombinations(4, :));
nkoutValues = length(koutValues);

nKoutSSEs = nan(nkoutValues, 1);
for i = 1:nkoutValues
    parameterIndices = parameterCombinations(4, :) == koutValues(i);
    
    minSSElocal = 0;
    for expIndex = 1:nExpResultsToCompare
        SSEs = squeeze(SSE_arr(:,expIndex,parameterIndices));
        minSSElocal = minSSElocal + min(min(SSEs));
    end
    nKoutSSEs(i) = minSSElocal;
end
figure();
loglog(koutValues, nKoutSSEs, 'k', 'LineWidth', 0.8);
hold on;
loglog([min(koutValues) max(koutValues)], [1 1] .* maxSSE, 'r--', 'LineWidth', 0.8);
xlabel('k_{out}');
ylabel('Total SSE');
xlim([min(koutValues) max(koutValues)]);
applyPaperFormatting

figFileName = [prefix '_profile_kout'];
if exportFigures
    print([folderName filesep 'pdf' filesep figFileName '.pdf'], '-dpdf');
    export_fig([folderName filesep 'png' filesep figFileName '.png'], '-q101', '-r300');
end

end

