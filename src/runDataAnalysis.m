% runDataAnalysis: Runs data analysis and plots profile figures, as well as
% figures related to the model.
%
%   Usage:
%   runDataAnalysis()

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

%% 1. Download & load dependencies
loadDependencies();

%% 2. Aggregate red/green profile CSV files into tables
aggregateData();

%% 3. Perform peak detection, compute lengths, and bin data
lengthAnalysis();

%% 4. Plot binned data
compareConditions();

%% 5. Simulate model & compute measurement model
% This is computationally intense and requires cluster usage.
% To run this, open and run src/modelFitting/runSampling.m
%%
addpath(['model']);
addpath(['model' filesep 'parameterEstimation']);

%% 6. Compare simulation results to data, and estimate parameters
compareToDataFastRestrictedAll();
compareToDataFastRestrictedWt();
compareToDataFastRestrictedbfa1bub2();
compareToDataFastRestrictedS63A();

%% 7. Export comparison plots
plotViolinPlotswtbfa1bub2();
plotViolinPlotswtS63A();
plotViolinPlotswtS63Abfa1bub2();