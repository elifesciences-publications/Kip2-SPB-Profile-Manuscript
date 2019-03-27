function aggregateData(fileName, binsToCompute, runSingleDebug, debugIndices)
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
indicesToRun = 1:optimizationWorkspace.nCombos;
if runSingleDebug
    indicesToRun = debugIndices;
end

parameterCombinations = optimizationWorkspace.parameterCombinations;

%% Compute profiles from MT data
jobs = [parallel.FevalFuture];
jobs = jobs(1:0, 1:0);

for i = indicesToRun
    p = parameterCombinations(:, i);
    if runSingleDebug
        aggregateResults(i, resultFolderCorrected, absResultFolder, binsToCompute, exportPlots, runSingleDebug, p);
    else
        jobs(end + 1) = parfeval(@aggregateResults, 0, i, resultFolderCorrected, absResultFolder,binsToCompute, exportPlots, runSingleDebug, p);
    end
end

if ~runSingleDebug
    h = waitbar(0, 'Computing...');
    nDone = 0;
    while nDone < optimizationWorkspace.nCombos
        try
            [idx] = fetchNext(jobs);
        catch ME
            errorIdx = 0;
            for idx = 1:optimizationWorkspace.nCombos
                if ~isempty(jobs(idx).Error)
                    errorIdx = idx;
                end
            end
            cancel([jobs(1:(errorIdx-1)) jobs((errorIdx+1):optimizationWorkspace.nCombos)]);
            warning('Job errored: %i', errorIdx);
            jobs(errorIdx).Error
            rethrow(ME)
        end
        nDone = nDone + 1;
        
        if mod(nDone, 20) == 0
            waitbar(nDone / optimizationWorkspace.nCombos, h, sprintf('Computing... %i / %i', nDone, optimizationWorkspace.nCombos));
        end
    end
    close(h) 
end
jobs = [];
%% Aggregate results
results = cell(1, optimizationWorkspace.nCombos);
for i = indicesToRun
    resultFile = load([absResultFolder.correctedMeanVm int2str(i) '.mat']);
    results{i} = resultFile.result;
end
save([folderPath filesep 'analyzed-' resultFolder extension], 'results', 'binsToCompute');
