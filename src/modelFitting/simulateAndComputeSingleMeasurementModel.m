function result = simulateAndComputeSingleMeasurementModel(absResultFolder, parameterIndex, currentParameters, tspan)

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

    if absResultFolder(end) ~= filesep
        absResultFolder = [absResultFolder filesep];
    end
    absVmResultFolder = [absResultFolder 'vm' filesep];

    currentFileName = '';
    result = simulateModel(currentParameters, parameterIndex, tspan, currentFileName);

    currentOutFileName = [absVmResultFolder int2str(parameterIndex) '.mat'];

    if ~isempty(currentParameters.binSizeLimit)
        result = computeMeasurementModel(result, currentOutFileName, currentParameters.binSizeLimit);
    else
        result = computeMeasurementModel(result, currentOutFileName);
    end

end

