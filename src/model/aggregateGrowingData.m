function aggregateGrowingData(fileName, binsToCompute)
% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

exportPlots = false;

optimizationWorkspace = load(fileName);
resultFolderCorrected = optimizationWorkspace.resultFolderCorrected;

[folderPath, resultFolder, extension] = fileparts(fileName);
resultFolderCorrectedMeanVm = [resultFolder '-corrected-mean-vm'];
resultFolderCorrectedMeanVmImages = [resultFolderCorrectedMeanVm filesep 'images'];

absPath = cd(cd([fileparts(mfilename('fullpath')) filesep '..' filesep '..' filesep 'simulationResults']));
toAbsPath = @(folder) [absPath filesep folder filesep];

absResultFolder = struct;
absResultFolder.correctedMeanVm       = toAbsPath(resultFolderCorrectedMeanVm);
absResultFolder.correctedMeanVmImages = toAbsPath(resultFolderCorrectedMeanVmImages);
%%
if exist(absResultFolder.correctedMeanVm, 'dir')
    movefile(absResultFolder.correctedMeanVm, [absResultFolder.correctedMeanVm(1:end-1) '-old-' datestr(clock, 30)]);
end
if ~exist(absResultFolder.correctedMeanVm, 'dir')
    mkdir(absResultFolder.correctedMeanVm(1:end-1));
end
mkdir(absResultFolder.correctedMeanVmImages(1:end-1));


%%
parameterCombinations = optimizationWorkspace.parameterCombinations;

%% Compute profiles from MT data
results = aggregateGrowingResults(resultFolderCorrected, absResultFolder, binsToCompute, exportPlots, parameterCombinations);

%% Aggregate results
save([folderPath filesep 'analyzed-' resultFolder extension], 'results', 'binsToCompute');
