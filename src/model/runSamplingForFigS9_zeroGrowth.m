%function resultFolder = runSampling(sampleMode) 
% Recommend running as script for easier debugging / resuming 

% runSampling: Script that samples the stochastic model and computes
% the measurement model in parameter space. Note that this procedure
% is fairly computationally intense, and was run on 400 cores to achieve
% fast compute times.

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)
clear

%%
startTime = datestr(clock, 30);
resultFolder = ['simResults-' startTime];

absPath = cd(cd([fileparts(mfilename('fullpath')) filesep '..' filesep '..' filesep 'simulationResults']));
toAbsPath = @(folder) [absPath filesep folder filesep];
toAbsPathNoFilesep = @(folder) [absPath filesep folder];

lengthsToCompute = (183-44):11:(349+44);

resultFolderCorrected = cell2table({NaN, resultFolder, toAbsPath(resultFolder), [toAbsPath(resultFolder) 'vm' filesep]}, 'VariableNames' , {'offset', 'folderName', 'absFolderName', 'absVmFolderName'});

%%
warning('off'); % Suppress warnings regarding table extension
for i = 1:length(lengthsToCompute)
    resultFolderCorrected.offset(i) = lengthsToCompute(i);
    resultFolderCorrected.folderName{i} = ['results' filesep resultFolder '-' int2str(lengthsToCompute(i))];
    resultFolderCorrected.absFolderName{i} = toAbsPath(resultFolderCorrected.folderName{i});
    resultFolderCorrected.absVmFolderName{i} = [resultFolderCorrected.absFolderName{i} 'vm' filesep]; 
end
warning('on'); % Reactivate

%% Create results folders
for i = 1:size(resultFolderCorrected,1)
    if ~exist(resultFolderCorrected.absFolderName{i}(1:end-1), 'dir')
        mkdir(resultFolderCorrected.absFolderName{i}(1:end-1));
    end
    if ~exist(resultFolderCorrected.absVmFolderName{i}(1:end-1), 'dir')
        mkdir(resultFolderCorrected.absVmFolderName{i}(1:end-1));
    end
end
%% Cluster profile
clusterProfile = 'SGE_BSSE'; % use 'local' if you want to run this on your own machine
nNodes = 400;

%% Debug options
runSingleDebug = true;    % Runs single-threaded, synchronously
runSingleParallel = false; % Runs single-threaded, asynchronously
debugIndices = 1;          % Which indices to run for debugging


%% Kip Concentrations
c_Kip2_wt = 61; % nM

volumeOfYeastCell = 50; % in um^3
nAvogadro         = 6.022e23; % 1/mol

% Assuming all protofilaments are symmetric, how many Kinesins do we have
% per protofilament?

nKip2Free = c_Kip2_wt * 1e-9; %               nmol/liter -> mol/liter
nKip2Free = nKip2Free * 1e-15; %              mol/liter  -> mol/um^3
nKip2Free = nKip2Free * nAvogadro; %          mol/um^3   -> 1/um^3
nKip2Free = nKip2Free * volumeOfYeastCell; %  1/um^3     -> 1
nKip2Free = nKip2Free / 13; %                 1 -> 1/protofilament
nKip2Free = round(nKip2Free) %                rounding


%% In vitro Kip2 lattice on rates (for reference)
k_on_per_um_MT_Kip2only       = 0.7; % per um MT, per nM Kip2, per minute, Roberts et al (2014), Figure 3 - figure supplement 
k_on_per_um_MT_Kip2_Bik1_Bim1 = 3.9; % per um MT, per nM Kip2, per minute, Roberts et al (2014), Figure 3 - figure supplement 2

robertsKonConverter = @(x) x / 13 / 1000 * 8 / 60; % per um MT, per nM Kip2, per min -> per site (PF, 8 nm sites), per second, per nM Kip2

k_on_Kip2only = robertsKonConverter(k_on_per_um_MT_Kip2only) % 1/(nM Kip2 * s) (per site)
k_on_Kip2_Bik1_Bim1 = robertsKonConverter(k_on_per_um_MT_Kip2_Bik1_Bim1) % 1/(nM Kip2 * s) (per site)

%% Stepping rate (Our data, fixed)
v_step_Kip2_vivo = 6.259 * 1000 / 60; % um/min -> nm/s, our data
k_step_Kip2_vivo = v_step_Kip2_vivo / 8; % nm/s -> sites/s


%% Off rate (Our data, fixed)
k_detach_Kip2_vivo = 0; % 1/s, de-facto in Xiuzhen's experimental data

%% Sampling time points
tend = 7200;
nPoints = 901;
tspan = linspace(3600, tend, nPoints);

%% Parameter space sampling settings
nKip2Free_range = [60:20:200 300:100:500]; % Kip2 molecules per protofilament
k_on_range      = logspace(-5, 0, 15); % per site (PF, 8 nm sites), per second, per nM Kip2
k_in_range      = logspace(-3, 2, 15); % per PF, per second, per nM Kip2
k_out_range     = linspace( 1, k_step_Kip2_vivo, 13); % per PF, per second

parameterCombinations = [80; 6.1054e-04; 0.3162; 3.0066];

%% Save parameter vector and setup for reproducibility
save([toAbsPathNoFilesep(resultFolder) '-start.mat']);

%%
nCombos = size(parameterCombinations, 2);
fprintf('\n# Parameter combinations to compute: %i\n', nCombos);

if ~runSingleDebug && ~runSingleParallel
    if isempty(gcp('nocreate')) && any(contains(parallel.clusterProfiles, clusterProfile))
        fprintf('\nFiring up %i cluster nodes at %s\n', nNodes,  datestr(clock, 30));
        parpool(clusterProfile, nNodes);
    end
end


%%
globalStart = tic;
fprintf('\nStarting simulations and measurement model computation at %s\n', datestr(clock, 30));

parameterIndices = 1:nCombos;
if runSingleParallel || runSingleDebug
    parameterIndices = debugIndices;
    
    warning(['Debug mode activated, only computing parameter indices ' mat2str(debugIndices)]);
    
    if runSingleParallel
        nCombos = length(debugIndices);
    end
end

jobs = [parallel.FevalFuture];
jobs = jobs(1:0, 1:0);

for parameterIndex = parameterIndices
    currentParams = parameterCombinations(:, parameterIndex);

    currentParameters                  = struct;
    
    currentParameters.lenghtsToCompute = lengthsToCompute;
    currentParameters.maxLength        = lengthsToCompute(1);
    currentParameters.nKip2Free        = currentParams(1);
    currentParameters.k_on_mt          = ones(1, currentParameters.maxLength) * currentParams(2);
    currentParameters.k_on_mt(1)       = currentParams(2) + currentParams(3);
    currentParameters.k_step_mt        = ones(1, currentParameters.maxLength) * k_step_Kip2_vivo;
    currentParameters.k_step_mt(end)   = 0;
    currentParameters.k_detach_mt      = zeros(1, currentParameters.maxLength) + k_detach_Kip2_vivo;
    currentParameters.k_detach_mt(end) = currentParams(4);
    currentParameters.reportRunLengths = false;

    if runSingleDebug
        runMeasurementModelForAllLengths(resultFolderCorrected, parameterIndex, currentParameters, tspan);
    else
        jobs(end + 1) = parfeval(@runMeasurementModelForAllLengths, 0, resultFolderCorrected, parameterIndex, currentParameters, tspan);
    end

end

if ~runSingleDebug
    h = waitbar(0, 'Computing...');
    nDone = 0;
    while nDone < nCombos
        [idx] = fetchNext(jobs);
        if ~isempty(jobs(idx).Error)
            idx
            jobs(idx).Error
            close(h) 
            cancel(jobs)
            error('Jobs errored - please investigate');
        end
        nDone = nDone + 1;
        
        if mod(nDone, 20) == 0
            waitbar(nDone / nCombos,h,sprintf('Computing... %i / %i', nDone, nCombos));
        end
    end
    close(h) 
end

globalRuntime.sim = toc(globalStart);
fprintf('\nSimulations and measurement model computation finished at %s\n', datestr(clock, 30));
%% Save simulation time
save([toAbsPathNoFilesep(resultFolder) '.mat']);

fprintf('\nStarting data binning at %s\n', datestr(clock, 30));
%% Bins to compute
binsToCompute = [
%    1.33 1.61
%    1.59 1.87
    1.85 2.14
%    2.13 2.41
%    2.39 2.67
%    2.65 2.94
];
%% Aggregate data in bins
aggregateData([toAbsPathNoFilesep(resultFolder) '.mat'], binsToCompute, runSingleDebug, debugIndices);
fprintf('\nData binning finished at %s\n', datestr(clock, 30));
%%
figureFolderAbsPath = cd(cd([fileparts(mfilename('fullpath')) filesep '..' filesep '..' filesep 'figures' filesep 'model']));
figureFolderAbsPathSep = [figureFolderAbsPath filesep];
export_fig([figureFolderAbsPathSep 'png' filesep 'FigS9_zeroGrowth.png'], '-q101', '-r300');
print([figureFolderAbsPathSep 'pdf' filesep 'FigS9_zeroGrowth.pdf'], '-dpdf');

%% Close parpool
if ~runSingleDebug
    delete(gcp('nocreate'));
end
