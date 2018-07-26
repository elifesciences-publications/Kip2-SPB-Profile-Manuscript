%function resultFolder = runSampling(sampleMode) 
% Recommend running as script for easier debugging / resuming 

% runSampling: Script that samples the stochastic model and computes
% the measurement model in parameter space. "sampleMode" switches
% between sampling modes used for Fig. 2B and 2C/S8. Two sets of parameter
% boundaries were used for Fig. 2C (see below). Note that this procedure
% is fairly computationally intense, and was run on 180-512 cores to achieve
% fast compute times.

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

clear
sampleMode = 'Fig2B';

switch sampleMode % Settings used for publication
    case 'Fig2C-1' % Used to generate "simResults-20180626T202346.mat"
        nNodes = 250;
        maxLength_range = 183; % sites (8 nm each) on PF
        nKip2Free_bounds = [100 1000];
        k_in_jitter_bounds = [0.5 1.5];
        k_out_jitter_bounds = [0.5 1.5];
        nSamples = 10000;
        
        bounds = [
            nKip2Free_bounds
            k_in_jitter_bounds
            k_out_jitter_bounds
        ];
        
        parameterCombinations(1, :) = maxLength_range * ones(1, nSamples);
        parameterCombinations([2 4 5], :) = repmat(bounds(:,1), [1, nSamples]) + lhsdesign(nSamples, 3)' .* repmat(bounds(:,2)-bounds(:,1), [1, nSamples]);
        parameterCombinations(3, :) = 4e-5 * ones(1, nSamples); % Enforce known on rate
        
        % Sample in approximately good region
        parameterCombinations(4, :) = parameterCombinations(4, :) .* 5.643 ./ (parameterCombinations(2, :) - 80.25);
        parameterCombinations(5, :) = parameterCombinations(5, :) .* (0.003 .* parameterCombinations(2, :) + 2.71);
        
        nPoints = 901; % Use 901 for sampling (4 Hz), 3601 for Fig. 2B (1 Hz)
    case 'Fig2C-2' % Used to generate "simResults-20180628T094544.mat"
        nNodes = 180;
        maxLength_range = 183; % sites (8 nm each) on PF
        nKip2Free_bounds = [81 600];
        k_in_jitter_bounds = [0.8 1.2];
        k_out_jitter_bounds = [0.8 1.2];
        nSamples = 2000;
        
        bounds = [
            nKip2Free_bounds
            k_in_jitter_bounds
            k_out_jitter_bounds
        ];
        
        parameterCombinations(1, :) = maxLength_range * ones(1, nSamples);
        parameterCombinations([2 4 5], :) = repmat(bounds(:,1), [1, nSamples]) + lhsdesign(nSamples, 3)' .* repmat(bounds(:,2)-bounds(:,1), [1, nSamples]);
        parameterCombinations(3, :) = 4e-5 * ones(1, nSamples); % Enforce known on rate
        
        % Sample in approximately good region
        parameterCombinations(4, :) = parameterCombinations(4, :) .* 5.643 ./ (parameterCombinations(2, :) - 80.25);
        parameterCombinations(5, :) = parameterCombinations(5, :) .* (0.003 .* parameterCombinations(2, :) + 2.71);
        
        nPoints = 901; % Use 901 for sampling (4 Hz), 3601 for Fig. 2B (1 Hz)
    case 'Fig2B' % Used to generate "simResults-20180706T153609.mat"
        nNodes = 4; % 6
        maxLength_range = [116 150 183 216 249 283]; % sites (8 nm each) on PF
        nKip2Free_range = 250; % molecules
        k_on_range      = 4e-5; % nM^-1 s^-1
        k_in_range      = 0.0358; % nM^-1 s^-1
        k_out_range     = 3.569; %3.731; % 1/s
        parameterCombinations = combvec(maxLength_range, nKip2Free_range, k_on_range, k_in_range, k_out_range);
        
        nPoints = 3601; % Use 901 for sampling (4 Hz), 3601 for Fig. 2B (1 Hz)
    case default
        error('Unknown mode. Supported modes: Fig2B, Fig2C-1, Fig2C-2');
end

%%
startTime = datestr(clock, 30);
resultFolder = ['simResults-' startTime];

absPath = cd(cd([fileparts(mfilename('fullpath')) filesep '..' filesep '..' filesep 'simulationResults']));
toAbsPath = @(folder) [absPath filesep folder filesep];
toAbsPathNoFilesep = @(folder) [absPath filesep folder];

offsetsToCompute = [-44 -33 -22 -11 0 11 22 33 44];
%%
warning('off'); % Suppress warnings regarding table extension
resultFolderCorrected = cell2table({NaN, resultFolder, toAbsPath(resultFolder), [toAbsPath(resultFolder) 'vm' filesep]}, 'VariableNames' , {'offset', 'folderName', 'absFolderName', 'absVmFolderName'});
resultFolderCorrected.offset(2) = NaN;
resultFolderCorrected.folderName{2} = [resultFolder '-binSizeCorrection'];
resultFolderCorrected.absFolderName{2} = toAbsPath(resultFolderCorrected.folderName{2});
resultFolderCorrected.absVmFolderName{2} = [resultFolderCorrected.absFolderName{2} 'vm' filesep]; 

for i = (1:length(offsetsToCompute))+2
    resultFolderCorrected.offset(i) = offsetsToCompute(i - 2);

    resultFolderCorrected.folderName{i} = [resultFolder '-corrected_' int2str(offsetsToCompute(i-2))];

    resultFolderCorrected.absFolderName{i} = toAbsPath(resultFolderCorrected.folderName{i});
    resultFolderCorrected.absVmFolderName{i} = [resultFolderCorrected.absFolderName{i} 'vm' filesep]; 
end
warning('on'); % Reactivate
for i = 1:size(resultFolderCorrected,1)
    mkdir(resultFolderCorrected.absFolderName{i}(1:end-1));
    mkdir(resultFolderCorrected.absVmFolderName{i}(1:end-1));
end

clusterProfile = 'local'; % use 'local' if you want to run this on your own machine

%% Debug options
runSingleDebug = false;    % Runs single-threaded, synchronously
runSingleParallel = false; % Runs single-threaded, asynchronously
debugIndices = 1;          % Which indices to run for debugging

%% Kip Concentrations
c_Kip2_wt = 30; % nM

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
% Result: nKip2Free = 69

%% In vitro Kip2 lattice on rate (fixed)
k_on_per_um_MT_Kip2_Bik1_Bim1 = 3.9; % per um MT, per nM Kip2, per minute, Roberts et al (2014), Figure 3 - figure supplement 2
robertsKonConverter = @(x) x / 13 / 1000 * 8 / 60; % per um MT, per nM Kip2, per min -> per site (PF, 8 nm sites), per second, per nM Kip2

k_on_Kip2_Bik1_Bim1 = robertsKonConverter(k_on_per_um_MT_Kip2_Bik1_Bim1); % 1/(nM Kip2 * s)

%% Stepping rate (Our data, fixed)
v_step_Kip2_vivo = 6.259 * 1000 / 60; % um/min -> nm/s, our data
k_step_Kip2_vivo = v_step_Kip2_vivo / 8; % nm/s -> sites/s

%% Off rate (fixed)
k_detach_Kip2 = 0; % 1/s, for justification see supplement

%% Set up samples
tend = 7200;
tspan = linspace(3600, tend, nPoints); % Use nPoints=901 for sampling (4 Hz), 3601 for Fig. 2B (1 Hz)

nCombos = size(parameterCombinations, 2);
fprintf('\n# Parameter combinations to compute: %i\n', nCombos);

%% Figure up cluster
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
end

jobs = [parallel.FevalFuture];
jobs = jobs(1:0, 1:0);

for parameterIndex = parameterIndices
    currentParams = parameterCombinations(:, parameterIndex);

    currentParameters                  = struct;
    currentParameters.maxLength        = currentParams(1); % MT protofilament length in sites (8 nm each)
    currentParameters.nKip2Free        = currentParams(2); % Number of Kip2 molecules available per protofilament
    currentParameters.k_on_mt          = ones(1, currentParameters.maxLength) * currentParams(3); % k_on (nM^-1 s^-1)
    currentParameters.k_on_mt(1)       = currentParams(3) + currentParams(4); % k_in + k_on (nM^-1 s^-1)
    currentParameters.k_step_mt        = ones(1, currentParameters.maxLength) * k_step_Kip2_vivo; % k_step (s^-1)
    currentParameters.k_step_mt(end)   = 0; % k_step at end (s^-1) - motor cannot step off
    currentParameters.k_detach_mt      = zeros(1, currentParameters.maxLength) + k_detach_Kip2; % k_off (s^-1)
    currentParameters.k_detach_mt(end) = currentParams(5); % k_out at end (s^-1)
    currentParameters.reportRunLengths = false; % optimize
    currentParameters.binSizeLimit     = currentParams(1)*0.008 + [-0.2 0.2]; % Peak must be +- pixel (plus half pixel where peak still is inside pixel)

    if runSingleDebug
        simulateAndComputeCorrectedMeasurementModel(resultFolderCorrected, parameterIndex, currentParameters, tspan);
    else
        jobs(end + 1) = parfeval(@simulateAndComputeCorrectedMeasurementModel, 0, resultFolderCorrected, parameterIndex, currentParameters, tspan);
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
%%
save([toAbsPathNoFilesep(resultFolder) '.mat']);
%%
fprintf('\n%s written\n', [toAbsPathNoFilesep(resultFolder) '.mat']);
%%
aggregateMeasurementModelResults([toAbsPathNoFilesep(resultFolder) '.mat'], resultFolder, runSingleDebug, debugIndices);
%%
if ~runSingleDebug
    delete(gcp('nocreate'));
end
