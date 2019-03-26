function result = runMeasurementModel(absResultFolder, parameterIndex, currentParameters, tspan)
    if absResultFolder(end) ~= filesep
        absResultFolder = [absResultFolder filesep];
    end
    absVmResultFolder = [absResultFolder 'vm' filesep];

    currentFileName = [absResultFolder int2str(parameterIndex) '.mat'];
    if ~exist(currentFileName, 'file')
        result = simulateModel(currentParameters, parameterIndex, tspan, currentFileName);
    end
    
    currentOutFileName = [absVmResultFolder int2str(parameterIndex) '.mat'];
    if ~exist(currentOutFileName, 'file')
        if isfield(currentParameters, 'binSizeLimit')
            result = computeMeasurementModel(result, currentOutFileName, currentParameters.binSizeLimit);
        else
            result = computeMeasurementModel(result, currentOutFileName);
        end
    end
end

