function result = runGrowingMeasurementModel(absResultFolder, currentParameters, tspan)
    if absResultFolder(end) ~= filesep
        absResultFolder = [absResultFolder filesep];
    end
    absVmResultFolder = [absResultFolder 'vm' filesep];

    result = cell(currentParameters.nReplicates, 1);
    parfor i = 1:currentParameters.nReplicates
        currentFileName = [absResultFolder int2str(i) '.mat'];
        if ~exist(currentFileName, 'file')
            result{i} = simulateGrowingModel(currentParameters, tspan, currentFileName);
        end
    end
    
    currentOutFileName = [absVmResultFolder 'simVm.mat'];
    if ~exist(currentOutFileName, 'file')
        result = computeGrowingMeasurementModel(result, currentParameters, currentOutFileName);
    end
end

