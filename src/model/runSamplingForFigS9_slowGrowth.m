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

%lengthsToCompute = (183-44):11:(349+44);
startLength = 200; %183-44;
endLength = 300; %349+44;

%startLength = 183-44;
%endLength = 349+44;

resultFolderCorrected = toAbsPath(resultFolder);

%% Create results folders
if ~exist(resultFolderCorrected, 'dir')
    mkdir(resultFolderCorrected);
end

if ~exist([resultFolderCorrected filesep 'vm'], 'dir')
    mkdir([resultFolderCorrected filesep 'vm']);
end

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


%% Stepping rate (Our data, fixed)
v_step_Kip2_vivo = 6.259 * 1000 / 60; % um/min -> nm/s, our data
k_step_Kip2_vivo = v_step_Kip2_vivo / 8; % nm/s -> sites/s


%% Off rate (Our data, fixed)
k_detach_Kip2_vivo = 0; % 1/s, de-facto in Xiuzhen's experimental data

%% Microtubule growth rate (Our data, fixed)
mtGrowthRate = 0.05; % 1/s (8 nm increments, each PF)
%mtGrowthRate = 2.875; % 1/s (8 nm increments, each PF)

%% Sampling time points
tend = 7200;
warmupTime = 3600;
nPoints = 901;
tspan = linspace(3600, tend, nPoints);

%% Parameter space sampling settings
nKip2Free = 80; % Kip2 molecules per protofilament
k_on      = 6.1054e-04; % per site (PF, 8 nm sites), per second, per nM Kip2
k_in      = 0.3162; % per PF, per second, per nM Kip2
k_out     = 3.0066; % per PF, per second

parameterCombinations = [nKip2Free; k_on; k_in; k_out];
nCombos = 1;

%% Save parameter vector and setup for reproducibility
save([toAbsPathNoFilesep(resultFolder) '-start.mat']);


%%
globalStart = tic;
fprintf('\nStarting simulations and measurement model computation at %s\n', datestr(clock, 30));


currentParams = parameterCombinations(:, 1);

currentParameters                  = struct;

currentParameters.nReplicates      = 100;
currentParameters.startLength      = startLength;
currentParameters.endLength        = endLength;
currentParameters.warmupTime       = warmupTime;
currentParameters.growthRate       = mtGrowthRate;
currentParameters.nKip2Free        = currentParams(1);
currentParameters.k_on_mt          = currentParams(2);
currentParameters.k_on_in_mt       = currentParams(2) + currentParams(3);
currentParameters.k_step_mt        = k_step_Kip2_vivo;
currentParameters.k_detach_mt      = 0;
currentParameters.k_detach_mt_end  = currentParams(4);
currentParameters.reportRunLengths = false;

runGrowingMeasurementModel(resultFolderCorrected, currentParameters, tspan);

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
aggregateGrowingData([toAbsPathNoFilesep(resultFolder) '.mat'], binsToCompute);
fprintf('\nData binning finished at %s\n', datestr(clock, 30));

%%
figureFolderAbsPath = cd(cd([fileparts(mfilename('fullpath')) filesep '..' filesep '..' filesep 'figures' filesep 'model']));
figureFolderAbsPathSep = [figureFolderAbsPath filesep];
export_fig([figureFolderAbsPathSep 'png' filesep 'FigS9_slowGrowth.png'], '-q101', '-r300');
print([figureFolderAbsPathSep 'pdf' filesep 'FigS9_slowGrowth.pdf'], '-dpdf');
