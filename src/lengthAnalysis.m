clear
close all
clc

%% Configure data folders
analyzedDataFolder = ['..' filesep 'analyzedData'];

dataFolder    = [analyzedDataFolder filesep 'aggregatedProfiles'];
dataFolderSep = [dataFolder filesep];

binnedDataFolder    = [analyzedDataFolder filesep 'aggregatedProfiles-binned'];
binnedDataFolderSep = [binnedDataFolder filesep];

if ~exist(binnedDataFolder, 'dir')
    mkdir(binnedDataFolder);
end

mergedDataFolder    = [analyzedDataFolder filesep 'aggregatedProfiles-merged'];
mergedDataFolderSep = [mergedDataFolder filesep];

if ~exist(mergedDataFolder, 'dir')
    mkdir(mergedDataFolder);
end

%% Configure figure folders
figureFolder    = ['..' filesep 'figures' filesep 'peakDetection'];
figureFolderSep = [figureFolder filesep];

figureFolderBinned  = ['..' filesep 'figures' filesep 'peakDetection-binned'];
figureFolderBinnedSep = [figureFolderBinned filesep];

%% Configure results folder
resultsFolder  = ['..' filesep 'analyzedData'];
resultsFolderSep = [resultsFolder filesep];

%% Scaling definitions
defaultScaling = struct;
defaultScaling.green.offset = 0;
defaultScaling.green.factor = 1;
defaultScaling.red.offset = 0;
defaultScaling.red.factor = 1;

%% Merge wt data (replicate strains, after inspecting data separately), while scaling red channel according to different exposure time
file1.green = [dataFolderSep '20180130_15100b_Kip2-3sfGFP_Spc42-mCherry-metaPhase-named-green.txt'];
file1.red   = [dataFolderSep '20180130_15100b_Kip2-3sfGFP_Spc42-mCherry-metaPhase-named-red.txt'];
scaling1 = defaultScaling;
file2.green = [dataFolderSep '20180130_15100_Kip2-3sfGFP_Spc42-mCherry-metaPhase-named-green.txt'];
file2.red   = [dataFolderSep '20180130_15100_Kip2-3sfGFP_Spc42-mCherry-metaPhase-named-red.txt'];
scaling2 = defaultScaling;
scaling2.red.factor = 2.5; % compensate for longer mCherry exposure
outFile.green = [mergedDataFolderSep '20180130_15100-merged-Kip2-3sfGFP-Spc42-mCherry-metaPhase-named-green.txt'];
outFile.red   = [mergedDataFolderSep '20180130_15100-merged-Kip2-3sfGFP-Spc42-mCherry-metaPhase-named-red.txt'];
mergeAndScaleDatasetsTwoColor(file1, scaling1, file2, scaling2, outFile)

%% Merge Kip3 distal data
kip3File1.green = [dataFolderSep '20180212_14590_Kip3-3sfGFP_Spc42-mCherry-Distal-named-green.txt'];
kip3File1.red   = [dataFolderSep '20180212_14590_Kip3-3sfGFP_Spc42-mCherry-Distal-named-red.txt'];
kip3File2.green = [dataFolderSep '20180212_14590b_Kip3-3sfGFP_Spc42-mCherry-Distal-named-green.txt'];
kip3File2.red   = [dataFolderSep '20180212_14590b_Kip3-3sfGFP_Spc42-mCherry-Distal-named-red.txt'];

outFileKip3.green = [mergedDataFolderSep '20180212_14590_merged_Kip3-3sfGFP_Spc42-mCherry-Distal-named-green.txt'];
outFileKip3.red   = [mergedDataFolderSep '20180212_14590_merged_Kip3-3sfGFP_Spc42-mCherry-Distal-named-red.txt'];
mergeAndScaleDatasetsTwoColor(kip3File1, defaultScaling, kip3File2, defaultScaling, outFileKip3)

%%
debug = false;
exportPlots = true;
exportResults = true;
opacity = 0.2;
ylims = [0.8 20]*1e4;
xlims = [0 6];
xlimsPlusEnd = xlims - 1;
xlimsSPB = [-6 1];
xlimsCentered = xlims - 1.5;

binBounds = [
    0.79 1.07
    1.05 1.34
    1.33 1.61
    1.59 1.87
    1.85 2.14
    2.13 2.41
    2.39 2.67
    2.65 2.94
    2.92 3.21
    3.19 3.47
];

binFileName = {
    '0_80-1_06'
    '1_06-1_33'
    '1_33-1_60'
    '1_60-1_86'
    '1_86-2_13'
    '2_13-2_40'
    '2_40-2_66'
    '2_66-2_93'
    '2_93-3_20'
    '3_20-3_46'
    };

binName = binFileName;
for i = 1:length(binName)
    binName{i} = strrep(binName{i},'_','.');
end

%%


dataSets = {};

dataSetTemplate = struct;
dataSetTemplate.greenThreshold = 1.47e4;
dataSetTemplate.redThreshold = 5e4;
dataSetTemplate.backgroundComputation = 'last';
dataSetTemplate.yAxisText = 'Kip2-3xsfGFP + Spc42-mCherry Fluorescence (AU)';
dataSetTemplate.exportBins = true;
% 
% dataSet = dataSetTemplate;
% dataSet.name = 'wt-20180130-15100b';
% dataSet.greenFile = file1.green;
% dataSet.redFile   = file1.red
% dataSets{end+1} = dataSet;
% 
% dataSet = dataSetTemplate;
% dataSet.name = 'wt-20180130-15100';
% dataSet.greenFile = file2.green;
% dataSet.redFile   = file2.green;
% dataSet.redScalingFactor = 2.5; % scale to account for different exposure time
% dataSet.redScalingOffset = 0;
% dataSets{end+1} = dataSet;
% 
dataSet = dataSetTemplate;
dataSet.name = 'wt-20180130-15100-merged';
dataSet.greenFile = outFile.green;
dataSet.redFile   = outFile.red;
dataSet.exportBins = true;
dataSets{end+1} = dataSet;

dataSet = dataSetTemplate;
dataSet.name = 'wt-20180130-15100b-Distal';
dataSet.greenFile = [dataFolderSep '20180130_15100b_Kip2-3sfGFP_Spc42-mCherry-Distal-named-green.txt']; 
dataSet.redFile   = [dataFolderSep '20180130_15100b_Kip2-3sfGFP_Spc42-mCherry-Distal-named-red.txt']; 
dataSet.exportBins = true;
dataSets{end+1} = dataSet;

dataSet = dataSetTemplate;
dataSet.name = 'wt-20180224-15100';
dataSet.greenFile = [dataFolderSep '20180224_15100_Kip2-3sfGFP_Spc42-mCherry-metaPhase-named-green.txt']; 
dataSet.redFile   = [dataFolderSep '20180224_15100_Kip2-3sfGFP_Spc42-mCherry-metaPhase-named-red.txt']; 
dataSet.exportBins = true;
dataSets{end+1} = dataSet;

dataSet = dataSetTemplate;
dataSet.name = 'wt-20180224-15100-Distal';
dataSet.greenFile = [dataFolderSep '20180224_15100_Kip2-3sfGFP_Spc42-mCherry-Distal-named-green.txt'];
dataSet.redFile   = [dataFolderSep '20180224_15100_Kip2-3sfGFP_Spc42-mCherry-Distal-named-red.txt']; 
dataSet.exportBins = true;
dataSets{end+1} = dataSet;

dataSet = dataSetTemplate;
dataSet.name = 'Kip2-S63A-20180224-15102';
dataSet.greenFile = [dataFolderSep '20180224_15102_Kip2-S63A-3sfGFP_Spc42-mCherry-metaPhase-named-green.txt']; 
dataSet.redFile   = [dataFolderSep '20180224_15102_Kip2-S63A-3sfGFP_Spc42-mCherry-metaPhase-named-red.txt']; 
dataSet.exportBins = true;
dataSets{end+1} = dataSet;

dataSet = dataSetTemplate;
dataSet.name = 'Kip2-S63A-20180224-15102-Distal';
dataSet.greenFile = [dataFolderSep '20180224_15102_Kip2-S63A-3sfGFP_Spc42-mCherry-Distal-named-green.txt']; 
dataSet.redFile   = [dataFolderSep '20180224_15102_Kip2-S63A-3sfGFP_Spc42-mCherry-Distal-named-red.txt']; 
dataSet.exportBins = true;
dataSets{end+1} = dataSet;

dataSet = dataSetTemplate;
dataSet.name = '20180212-14590b-Kip3';
dataSet.greenFile = [dataFolderSep '20180212_14590b_Kip3-3sfGFP_Spc42-mCherry-metaPhase-named-green.txt']; 
dataSet.redFile   = [dataFolderSep '20180212_14590b_Kip3-3sfGFP_Spc42-mCherry-metaPhase-named-red.txt'];
dataSet.exportBins = true;
dataSet.yAxisText = 'Kip3-3xsfGFP + Spc42-mCherry Fluorescence (AU)';
dataSet.backgroundComputation = 'first'; % Need to do background correction past plus end (rather than SPB) due to fluorescence in the nucleus
dataSets{end+1} = dataSet;

dataSet = dataSetTemplate;
dataSet.name = '20180212-14590-Kip3-Distal-all';
dataSet.greenFile = outFileKip3.green; 
dataSet.redFile   = outFileKip3.red;
dataSet.exportBins = true;
dataSet.yAxisText = 'Kip3-3xsfGFP + Spc42-mCherry Fluorescence (AU)';
dataSet.backgroundComputation = 'first'; % Need to do background correction past plus end (rather than SPB) due to fluorescence in the nucleus
dataSets{end+1} = dataSet;

%% Create binned datasets
for i = 1:length(dataSets)
    if isfield(dataSets{i}, 'exportBins') && dataSets{i}.exportBins
        for j = 1:size(binBounds, 1)
            dataSet = dataSets{i}; % Inherit dataset properties
            dataSet.name = [dataSets{i}.name '-bin-' binName{j}];
            dataSet.greenFile = [binnedDataFolderSep dataSets{i}.name '-bin_' binFileName{j} '-green.txt'];
            dataSet.redFile   = [binnedDataFolderSep dataSets{i}.name '-bin_' binFileName{j} '-red.txt'];
            dataSet.doFlatness = true;
            dataSets{end+1} = dataSet;    
        end
    end
end
%%
nConditions = length(dataSets);

conditionResults = {};
for conditionIndex = 1:nConditions
    currentDataSet = dataSets{conditionIndex};
    close all
    
    greenDataTable = readtable(currentDataSet.greenFile);
    redDataTable   = readtable(currentDataSet.redFile);
    
    lengthVector = greenDataTable{:,1}; % in micrometers
    
    keepColumn = [];
    for i = 2:size(greenDataTable, 2)
        if ~iscell(greenDataTable{1, i})
            if ~isfield(currentDataSet, 'omitColumns') || isempty(intersect(currentDataSet.omitColumns,i))
                keepColumn = [keepColumn i];
            end
        end
    end
    
    % Detect duplicates, if any, and remove them from the analysis
    greenIntensities = greenDataTable{:,keepColumn}; % in AU
    redIntensities   = redDataTable{:,keepColumn};
    
    intensitiesForUniqueComparison = greenIntensities;
    intensitiesForUniqueComparison(isnan(intensitiesForUniqueComparison)) = -Inf;
    [intensitiesForUniqueComparison, keptCols, ~] = unique(intensitiesForUniqueComparison', 'stable', 'rows');
    intensitiesForUniqueComparison = intensitiesForUniqueComparison';
    intensitiesForUniqueComparison(isinf(intensitiesForUniqueComparison)) = NaN;
    keepColumn = keepColumn(keptCols);
    if size(greenIntensities, 2) ~= size(intensitiesForUniqueComparison, 2)
        warning(sprintf('%i duplicate profiles detected! Disregarding duplicates.', size(greenIntensities, 2) - size(intensitiesForUniqueComparison, 2)));
        greenIntensities = intensitiesForUniqueComparison;
        redIntensities = redIntensities(:, keptCols);
    end
    
    if isfield(currentDataSet, 'greenScalingFactor')
        greenIntensities = greenIntensities .* currentDataSet.greenScalingFactor + currentDataSet.greenScalingOffset;
    end
    
    if isfield(currentDataSet, 'redScalingFactor')
        redIntensities = redIntensities .* currentDataSet.redScalingFactor + currentDataSet.redScalingOffset;
    end
    
    mtPresent = ~isnan(greenIntensities);
    isMT = ~all(~mtPresent); % Exclude empty columns
    nMicrotubules = sum(isMT);
    
    if nMicrotubules < 2
        warning('Not enough data for current dataset, skipping');
        continue;
    end
    
    greenIntensities = greenIntensities(:, isMT);
    redIntensities   = redIntensities(:, isMT);
    mtIndices = find(isMT);
    keepColumn = keepColumn(isMT);
    
    %% Find peaks
    greenResults = {};
    redResults = {};
    first = [];
    last = [];
    
    greenThreshold = currentDataSet.greenThreshold;
    redThreshold   = currentDataSet.redThreshold;

    if debug
        figure(1);
    end
    for i = 1:nMicrotubules
        if debug
            clf
            set(gca, 'YScale', 'log');
            hold on;
        end
        thisMTGreenIntensity = greenIntensities(:, i);
        thisMTRedIntensity   = redIntensities(:, i);
        thisMTlength = find(isnan(thisMTGreenIntensity), 1) - 1; % Last MT index
        if find(isnan(thisMTRedIntensity), 1) - 1 ~= thisMTlength
            error('Red and green MT length must be the same!');
        end
        if isempty(thisMTlength)
            thisMTlength = length(thisMTGreenIntensity);
        end
        
        if debug
            plot(lengthVector, thisMTGreenIntensity, 'g');
            
            plot(lengthVector, thisMTRedIntensity, 'r');
            
        end

        greenResult = struct;
        [greenResult.pks, greenResult.locs, greenResult.w, greenResult.p] = findpeaks(thisMTGreenIntensity);
        
        redResult = struct;
        [redResult.pks, redResult.locs, redResult.w, redResult.p] = findpeaks(thisMTRedIntensity);

        firstIndex = 1;
        lastIndex = length(redResult.pks);
        
        if debug
            plot(lengthVector(greenResult.locs(firstIndex)), greenResult.pks(firstIndex), '*');
            plot(lengthVector(redResult.locs(lastIndex)), redResult.pks(lastIndex), 'x');
            plot(lengthVector([1 thisMTlength]), greenThreshold*[1 1], 'g--');
            plot(lengthVector([1 thisMTlength]), redThreshold*[1 1], 'r--');
            title([int2str(i) ': ' greenDataTable.Properties.VariableNames{mtIndices(i)+1}])
        end

        if greenResult.pks(firstIndex) < greenThreshold 
            warning('First peak is low - please inspect data for correctness. Proceeding to next peak that is over %g!', greenThreshold);
            firstIndex = find(greenResult.pks >= greenThreshold, 1);
            if isempty(firstIndex) || greenResult.locs(firstIndex) > 0.5 * thisMTlength
                firstIndex = 1;
                warning('Next peak not within first half of profile, using first peak even though it is low.');
            end
            if debug
                plot(lengthVector(greenResult.locs(firstIndex)), greenResult.pks(firstIndex), 'or');
            end
        end

        if redResult.pks(lastIndex) < redThreshold
            warning('Last peak is low - please inspect data for correctness. Proceeding to next peak that is over %g!', redThreshold);
            lastIndex = find(redResult.pks >= redThreshold, 1, 'last');
            if isempty(lastIndex) || redResult.locs(lastIndex) < 0.5 * thisMTlength
                lastIndex = length(redResult.pks);
                 warning('Next peak not within last half of profile, using last peak even though it is low.');
            end
            if debug
                plot(lengthVector(redResult.locs(lastIndex)), redResult.pks(lastIndex), 'or');
            end
        end
        
        if greenResult.locs(firstIndex) == redResult.locs(lastIndex)
            i+1
            warning('Microtubule of length 0 - skipping!');
        end
        if debug
            pause;
        end
        first = [first; greenResult.pks(firstIndex) greenResult.locs(firstIndex) greenResult.w(firstIndex) greenResult.p(firstIndex)];
        last = [last; redResult.pks(lastIndex) redResult.locs(lastIndex) redResult.w(lastIndex) redResult.p(lastIndex)];
        greenResults{end + 1} = greenResult;
        redResults{end + 1} = redResult;
       
    end
    
    %%
    figure(2);
    lengths = last(:, 2) - first(:, 2);
    lengthsUM = lengths .* 4/30;
    histogram(lengthsUM);
    
    if debug
        figure(3);
        subplot(1,4,1);
        histogram(first(:,1));
        subplot(1,4,2);
        histogram(first(:,2));
        subplot(1,4,3);
        histogram(first(:,3));
        subplot(1,4,4);
        histogram(first(:,4));
    end
    
    %% Compute mean & median with 95% CI
    fprintf('\n%s\n------------------------------\n', currentDataSet.name);
    totalLengths = sum(mtPresent(:, isMT));
    totalLengthsUM = totalLengths .* 4/30;
    meanLength = mean(totalLengthsUM);
    meanLengthCI = bootci(10000,@mean,totalLengthsUM);
    fprintf('Mean profile length: \t%g um \t[%g - %g]\n', meanLength, meanLengthCI(1), meanLengthCI(2));
    medianLength = median(totalLengthsUM);
    medianLengthCI = bootci(10000,@median,totalLengthsUM);
    fprintf('Median profile length: \t%g um \t[%g - %g]\n', medianLength, medianLengthCI(1), medianLengthCI(2));
    meanLength = mean(lengthsUM(lengthsUM > 0));
    meanLengthCI = bootci(10000,@mean,lengthsUM(lengthsUM > 0));
    fprintf('Mean peak-to-peak length: \t%g um \t[%g - %g]\n', meanLength, meanLengthCI(1), meanLengthCI(2));
    medianLength = median(lengthsUM(lengthsUM > 0));
    medianLengthCI = bootci(10000,@median,lengthsUM(lengthsUM > 0));
    fprintf('Median peak-to-peak length: \t%g um \t[%g - %g]\n', medianLength, medianLengthCI(1), medianLengthCI(2));
    
    
    %%
    offset = max(first(:,2));
    greenIntensitiesShifted = nan(size(greenIntensities, 1) + offset, nMicrotubules);
    redIntensitiesShifted   = nan(size(redIntensities  , 1) + offset, nMicrotubules);
    
    lengthVectorShifted = ((-(offset-1):length(lengthVector)) * (4/30))';
    
    for i = 1:nMicrotubules
        greenIntensitiesShifted((offset - first(i,2) + 1):(size(greenIntensities,1) + offset - first(i,2)), i) = greenIntensities(:, i);
        redIntensitiesShifted((offset - first(i,2) + 1):(size(redIntensities,1) + offset - first(i,2)), i) = redIntensities(:, i);
    end
    
    offsetSPB = max(last(:,2));
    greenIntensitiesSPBShifted = nan(size(greenIntensities, 1) + offsetSPB, nMicrotubules);
    redIntensitiesSPBShifted   = nan(size(redIntensities,   1) + offsetSPB, nMicrotubules);
    lengthVectorSPBShifted = ((-(offsetSPB-1):length(lengthVector)) * (4/30))';
    
    for i = 1:nMicrotubules
        greenIntensitiesSPBShifted((offsetSPB - last(i,2) + 1):(size(greenIntensities,1) + offsetSPB - last(i,2)), i) = greenIntensities(:, i);
        redIntensitiesSPBShifted((offsetSPB - last(i,2) + 1):(size(redIntensities,1) + offsetSPB - last(i,2)), i) = redIntensities(:, i);
    end
    
    meanGreenIntensity = mean(greenIntensities,2, 'omitnan');
    meanRedIntensity = mean(redIntensities,2,'omitnan');
    
    nMTs = sum(~isnan(greenIntensities),2);
    
    semGreenIntensity  = std(greenIntensities,0,2, 'omitnan')./sqrt(nMTs);
    semGreenIntensity(nMTs < 2) = NaN;
    meanGreenIntensity(nMTs < 2) = NaN;
    
    semRedIntensity  = std(redIntensities,0,2, 'omitnan')./sqrt(nMTs);
    semRedIntensity(nMTs < 2) = NaN;
    meanRedIntensity(nMTs < 2) = NaN;    
    
    meanGreenIntensityShifted = mean(greenIntensitiesShifted,2, 'omitnan');
    nMTsShifted = sum(~isnan(greenIntensitiesShifted),2);
    semGreenIntensityShifted = std(greenIntensitiesShifted,0,2, 'omitnan')./sqrt(nMTsShifted);
    semGreenIntensityShifted(nMTsShifted < 2) = NaN;
    meanGreenIntensityShifted(nMTsShifted < 2) = NaN;

    meanRedIntensityShifted = mean(redIntensitiesShifted,2, 'omitnan');
    semRedIntensityShifted = std(redIntensitiesShifted,0,2, 'omitnan')./sqrt(nMTsShifted);
    semRedIntensityShifted(nMTsShifted < 2) = NaN;
    meanRedIntensityShifted(nMTsShifted < 2) = NaN;
    
    
    meanGreenIntensitySPBShifted = mean(greenIntensitiesSPBShifted,2, 'omitnan');
    nMTsSPBShifted = sum(~isnan(greenIntensitiesSPBShifted),2);
    stdGreenIntensitySPBShifted = std(greenIntensitiesSPBShifted,0,2, 'omitnan');
    semGreenIntensitySPBShifted = stdGreenIntensitySPBShifted./sqrt(nMTsSPBShifted);
    semGreenIntensitySPBShifted(nMTsSPBShifted < 2) = NaN;
    meanGreenIntensitySPBShifted(nMTsSPBShifted < 2) = NaN;
    
    meanRedIntensitySPBShifted = mean(redIntensitiesSPBShifted,2, 'omitnan');
    stdRedIntensitySPBShifted = std(redIntensitiesSPBShifted,0,2, 'omitnan');
    semRedIntensitySPBShifted = stdRedIntensitySPBShifted./sqrt(nMTsSPBShifted);
    semRedIntensitySPBShifted(nMTsSPBShifted < 2) = NaN;
    meanRedIntensitySPBShifted(nMTsSPBShifted < 2) = NaN;
    
    maxMeanOffset = max(first(:,2) + last(:,2));
    minMeanOffset = min(first(:,2) + last(:,2));
    meanMeanOffset = round(mean(last(:,2) - first(:,2)));
    
    greenIntensitiesCenterShifted = nan(2*size(greenIntensities,1) - 1 + maxMeanOffset, nMicrotubules);
    redIntensitiesCenterShifted = nan(2*size(redIntensities,1) - 1 + maxMeanOffset, nMicrotubules);
    lengthVectorCenterShifted = (((-(maxMeanOffset-1):(2*size(greenIntensities,1) - 1))+meanMeanOffset) * (2/30))'; 
    
    for i = 1:nMicrotubules
        greenIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 1):2:(size(greenIntensities,1)*2-1 + maxMeanOffset - (first(i,2) + last(i,2))), i) = greenIntensities(:, i);
        greenIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 2):2:(size(greenIntensities,1)*2-1 + maxMeanOffset - (first(i,2) + last(i,2))), i) = ...
            (greenIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 1):2:(size(greenIntensities,1)*2-2 + maxMeanOffset - (first(i,2) + last(i,2))), i) + ...
            greenIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 3):2:(size(greenIntensities,1)*2-1 + maxMeanOffset - (first(i,2) + last(i,2))), i)) * 0.5;
        
        redIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 1):2:(size(redIntensities,1)*2-1 + maxMeanOffset - (first(i,2) + last(i,2))), i) = redIntensities(:, i);
        redIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 2):2:(size(redIntensities,1)*2-1 + maxMeanOffset - (first(i,2) + last(i,2))), i) = ...
            (redIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 1):2:(size(redIntensities,1)*2-2 + maxMeanOffset - (first(i,2) + last(i,2))), i) + ...
            redIntensitiesCenterShifted((maxMeanOffset - (first(i,2) + last(i,2)) + 3):2:(size(redIntensities,1)*2-1 + maxMeanOffset - (first(i,2) + last(i,2))), i)) * 0.5;
    end
    
    meanGreenIntensityCenterShifted = mean(greenIntensitiesCenterShifted,2, 'omitnan');
    stdGreenIntensityCenterShifted = std(greenIntensitiesCenterShifted,0,2, 'omitnan');
    nMTsCenterShifted = sum(~isnan(greenIntensitiesCenterShifted),2);
    semGreenIntensityCenterShifted = stdGreenIntensityCenterShifted./sqrt(nMTsCenterShifted);
    
    semGreenIntensityCenterShifted(nMTsCenterShifted < 2) = NaN;
    meanGreenIntensityCenterShifted(nMTsCenterShifted < 2) = NaN;
    
    
    meanRedIntensityCenterShifted = mean(redIntensitiesCenterShifted,2, 'omitnan');
    stdRedIntensityCenterShifted = std(redIntensitiesCenterShifted,0,2, 'omitnan');
    semRedIntensityCenterShifted = stdRedIntensityCenterShifted./sqrt(nMTsCenterShifted);
    
    semRedIntensityCenterShifted(nMTsCenterShifted < 2) = NaN;
    meanRedIntensityCenterShifted(nMTsCenterShifted < 2) = NaN;
    
    %%
    figure(8);
    clf;
    subplot(1,4,1);
    hold on;
    
    uniqueLengths = unique(lengths(lengths > 0));
    nUniqueLengths = length(uniqueLengths);
    
    for i=1:nMicrotubules
        if lengths(i) == 0
            plot(lengthVector,greenIntensities(:, i), 'y');
            plot(lengthVector,redIntensities(:, i), 'y');
        else
            foo = plot(lengthVector,greenIntensities(:, i), 'g');
            foo.Color(4) = opacity;
            
            foo = plot(lengthVector,redIntensities(:, i), 'r');
            foo.Color(4) = opacity;
        end
    end
    
    shadedErrorBar(lengthVector,meanGreenIntensity, semGreenIntensity);
    shadedErrorBar(lengthVector,meanRedIntensity, semRedIntensity);
    xlabel('Distance (raw, {\mu}m)');
    ylabel(currentDataSet.yAxisText);
    %xlim([0 lengthVector(find(~all(isnan(intensities),2),1,'last'))]);
    xlim(xlims);
    ylim(ylims);
    set(gca, 'YScale', 'log');
    grid on;
    
    subplot(1,4,2);
    hold on;
    
    for i = 1:nMicrotubules
        if lengths(i) == 0
            plot(lengthVectorShifted,greenIntensitiesShifted(:, i), 'y');
            plot(lengthVectorShifted,redIntensitiesShifted(:, i), 'y');
        else
            foo = plot(lengthVectorShifted, greenIntensitiesShifted(:, i), 'g');
            foo.Color(4) = opacity;
            
            foo = plot(lengthVectorShifted, redIntensitiesShifted(:, i), 'r');
            foo.Color(4) = opacity;
        end
        
    end
%     fooShifted = plot(lengthVectorShifted,intensitiesShifted, 'b');
%     for i=1:nMicrotubules
%         fooShifted(i).Color(4) = opacity;
%     end
    shadedErrorBar(lengthVectorShifted, meanGreenIntensityShifted, semGreenIntensityShifted);
    shadedErrorBar(lengthVectorShifted, meanRedIntensityShifted, semRedIntensityShifted);
    %xlim([lengthVectorShifted(find(sum(~isnan(intensitiesShifted),2) >= 3,1)) lengthVectorShifted(find(sum(~isnan(intensitiesShifted),2) >= 3,1,'last'))]);
    xlim(xlimsPlusEnd);
    ylim(ylims);
    xlabel('Distance from plus end ({\mu}m)');
    applyPaperFormatting
    set(gca, 'YScale', 'log');
    grid on;


    subplot(1,4,3);
    hold on;

    for i = 1:nMicrotubules
        if lengths(i) == 0
            plot(lengthVectorSPBShifted, greenIntensitiesSPBShifted(:, i), 'y');
            plot(lengthVectorSPBShifted, redIntensitiesSPBShifted(:, i), 'y');
        else
            foo = plot(lengthVectorSPBShifted, greenIntensitiesSPBShifted(:, i), 'g');
            foo.Color(4) = opacity;
            
            foo = plot(lengthVectorSPBShifted, redIntensitiesSPBShifted(:, i), 'r');
            foo.Color(4) = opacity;            
        end
        
    end

    hold on;
    shadedErrorBar(lengthVectorSPBShifted, meanGreenIntensitySPBShifted, semGreenIntensitySPBShifted);
    shadedErrorBar(lengthVectorSPBShifted, meanRedIntensitySPBShifted, semRedIntensitySPBShifted);
    
    xlim(xlimsSPB);
    ylim(ylims);
    xlabel('Distance from SPB ({\mu}{m})');
    applyPaperFormatting
    set(gca, 'YScale', 'log');
    grid on;
    
    subplot(1,4,4);
    hold on;
    for i = 1:nMicrotubules
        if lengths(i) == 0
            plot(lengthVectorCenterShifted, greenIntensitiesCenterShifted(:, i), 'y');
            plot(lengthVectorCenterShifted, redIntensitiesCenterShifted(:, i), 'y');
        else
            foo = plot(lengthVectorCenterShifted, greenIntensitiesCenterShifted(:, i), 'g');
            foo.Color(4) = opacity;
            
            foo = plot(lengthVectorCenterShifted, redIntensitiesCenterShifted(:, i), 'r');
            foo.Color(4) = opacity;
        end
    end
    shadedErrorBar(lengthVectorCenterShifted, meanGreenIntensityCenterShifted, semGreenIntensityCenterShifted);% s[upperQ lowerQ]);%stdIntensitySPBShifted);
    shadedErrorBar(lengthVectorCenterShifted, meanRedIntensityCenterShifted, semRedIntensityCenterShifted);% s[upperQ lowerQ]);%stdIntensitySPBShifted);
    xlim(xlimsCentered);
    ylim(ylims);
    xlabel('Distance (center-aligned, {\mu}m)');
    applyPaperFormatting
    set(gca, 'YScale', 'log');
    grid on;

    foo = gcf;
    

    foo.Color = 'white';
    foo.Position(3) = 2*foo.Position(3);
    movegui(foo, 'center')

    if exportPlots
        export_fig([figureFolderSep 'png' filesep currentDataSet.name '-all.png'], '-q101', '-r300');
        fig.PaperPositionMode = 'manual';
        orient(fig,'landscape')
        print([figureFolderSep 'pdf' filesep currentDataSet.name '-all.pdf'], '-dpdf');
    end
    
    % Background intensity extraction 
    switch currentDataSet.backgroundComputation
        case 'last'
            % Last pixel past the SPB is background
            backgroundIntensities = nan(nMicrotubules, 1);
            for i = 1:nMicrotubules
                backgroundIntensities(i) = greenIntensities(totalLengths(i),i);
            end
        case 'first'
            % First pixel before plus end is background
            backgroundIntensities = nan(nMicrotubules, 1);
            for i = 1:nMicrotubules
                backgroundIntensities(i) = greenIntensities(1,i);
            end
        case 'both'
            % Both of the above are background
            backgroundIntensities = nan(2*nMicrotubules, 1);
            for i = 1:nMicrotubules
                backgroundIntensities(2*i-1) = greenIntensities(1,i);
                backgroundIntensities(2*i) = greenIntensities(totalLengths(i),i);
            end            
    end
    
    if isfield(currentDataSet, 'doFlatness') && currentDataSet.doFlatness
        figure(22);
        hold on;
        for i = 1:nMicrotubules
            if lengths(i) == 0
                plot(lengthVectorCenterShifted, greenIntensitiesCenterShifted(:, i), 'y');
            else
                foo = plot(lengthVectorCenterShifted, greenIntensitiesCenterShifted(:, i), 'g');
                foo.Color(4) = opacity;
                
                foo = plot(lengthVectorCenterShifted, redIntensitiesCenterShifted(:, i), 'r');
                foo.Color(4) = opacity;
            end
            

        end
        H = shadedErrorBar(lengthVectorCenterShifted, meanGreenIntensityCenterShifted, semGreenIntensityCenterShifted);
        H = shadedErrorBar(lengthVectorCenterShifted, meanRedIntensityCenterShifted, semRedIntensityCenterShifted);
        
        xlim(xlimsCentered);
        ylim(ylims);
        xlabel('Distance from plus end ({\mu}m)');
        ylabel(currentDataSet.yAxisText);
        applyPaperFormatting
        set(gca, 'YScale', 'log');
        grid on;
        
        H_BG  = shadedErrorBar(xlimsCentered, [1 1]*mean(backgroundIntensities), [1 1]*std(backgroundIntensities)./sqrt(length(backgroundIntensities)), 'lineProps', {'--', 'Color', [0.5 0.5 0.5]});
        H_SPB = plot([meanMeanOffset-0.5 meanMeanOffset-0.5]*(4/30), ylims,  '--r');


        flatnessResult = struct;
        legend([H.mainLine, H.edge(1), H_BG.mainLine], {['Mean, n = ' int2str(nMicrotubules)], 'SEM', 'Mean Background'}, 'Location', 'NorthWest');

        legend boxoff;

        foo = gcf;
        foo.Color = 'white';
        
        movegui(foo, 'center')
        if exportPlots
            export_fig([figureFolderBinnedSep 'png' filesep currentDataSet.name '-centerAligned.png'], '-q101', '-r300');
            print([figureFolderBinnedSep 'pdf' filesep currentDataSet.name '-centerAligned.pdf'], '-dpdf');
        end
    end

    %%

    
    condResult = struct;
    condResult.dataSet = currentDataSet;
    condResult.totalLengths = totalLengths;
    condResult.totalLengthsUM = totalLengthsUM;
    condResult.lengthsUM = lengthsUM;
    condResult.lengths = lengths;
    condResult.condition = currentDataSet.name;
    condResult.lengthVector = lengthVector;
    condResult.greenIntensities = greenIntensities;
    condResult.redIntensities   = redIntensities;
    condResult.offset = offset;
    condResult.offsetSPB = offsetSPB;
    condResult.greenIntensitiesShifted = greenIntensitiesShifted;
    condResult.redIntensitiesShifted = redIntensitiesShifted;
    condResult.lengthVectorShifted = lengthVectorShifted;
    condResult.greenIntensitiesSPBShifted = greenIntensitiesSPBShifted;
    condResult.redIntensitiesSPBShifted = redIntensitiesSPBShifted;
    condResult.lengthVectorSPBShifted = lengthVectorSPBShifted;
    condResult.backgroundIntensities = backgroundIntensities;
    
    condResult.lengthVectorCenterShifted = lengthVectorCenterShifted;
    condResult.greenIntensitiesCenterShifted = greenIntensitiesCenterShifted;
    condResult.redIntensitiesCenterShifted = redIntensitiesCenterShifted;
    condResult.meanMeanOffset = meanMeanOffset;
    
    condResult.greenPeaks = greenResults;
    condResult.redPeaks   = redResults;
    condResult.cellNames = greenDataTable.Properties.VariableNames(keepColumn);
    
    if isfield(currentDataSet, 'doFlatness') && currentDataSet.doFlatness
        condResult.flatnessResult = flatnessResult;
    end
    conditionResults{end + 1} = condResult;
    
    
    if isfield(currentDataSet, 'exportBins') && currentDataSet.exportBins
        prefix = [binnedDataFolderSep currentDataSet.name '-'];

        for i = 1:size(binBounds, 1)
            indices = find(condResult.lengthsUM > binBounds(i, 1) & condResult.lengthsUM < binBounds(i, 2));
            binnedCondition = createTwoColorSubset(condResult, indices, [condResult.condition '-bin_' binName{i}]);
            exportTwoColorCondtionResultToData(binnedCondition, [prefix 'bin_' binFileName{i}]);   
        end
    end
    
    if ~exportResults
        pause
    end
end



%%
if exportResults
    save([resultsFolderSep 'binnedProfiles.mat'], 'conditionResults');
end
