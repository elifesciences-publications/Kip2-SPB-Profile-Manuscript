function simulateAndComputeCorrectedMeasurementModel(resultFolderCorrected, parameterIndex, currentParameters, tspan)
    currentParametersInitial = currentParameters;
    currentParametersInitial.binSizeLimit = [];
    result = simulateAndComputeSingleMeasurementModel(resultFolderCorrected.absFolderName{1}, parameterIndex, currentParametersInitial, tspan);
    
    correctionOffset = round(result.model.meanPlusEndLengthOffset / 0.008);
    
    if correctionOffset > 150
        warning(sprintf('Correction offset capped at 150 instead of %i', correctionOffset));
        correctionOffset = 150;
    elseif correctionOffset < -currentParameters.maxLength/2
        warning(sprintf('Correction offset capped at %i instead of %i', -round(currentParameters.maxLength/2), correctionOffset));
        correctionOffset = -round(currentParameters.maxLength/2);
    end
    
    
    currentParametersInitialBinSizeLimit = updateParameters(currentParameters, correctionOffset);
    resultBinSizeLimit = simulateAndComputeSingleMeasurementModel(resultFolderCorrected.absFolderName{2}, parameterIndex, currentParametersInitialBinSizeLimit, tspan);
    
    correctionOffsetBinSizeLimit = round(resultBinSizeLimit.model.meanPlusEndLengthOffset / 0.008);
    
    if correctionOffsetBinSizeLimit > 150
        warning(sprintf('Correction offset (2) capped at 150 instead of %i', correctionOffsetBinSizeLimit));
        correctionOffsetBinSizeLimit = 150;
    elseif correctionOffsetBinSizeLimit < -currentParameters.maxLength/2
        warning(sprintf('Correction offset (2) capped at %i instead of %i', -round(currentParameters.maxLength/2), correctionOffsetBinSizeLimit));
        correctionOffsetBinSizeLimit = -round(currentParameters.maxLength/2);
    end
    
    function currentParametersUpdated = updateParameters(currentParameters, offset)
        currentParametersUpdated = currentParameters;
        currentParametersUpdated.maxLength   = currentParameters.maxLength + offset;
        
        if offset > 0
            currentParametersUpdated.k_on_mt     = [currentParameters.k_on_mt currentParameters.k_on_mt(end) * ones(1, offset)];
            currentParametersUpdated.k_step_mt   = [currentParameters.k_step_mt(1) * ones(1, offset) currentParameters.k_step_mt];
            currentParametersUpdated.k_detach_mt = [currentParameters.k_detach_mt(1) * ones(1, offset) currentParameters.k_detach_mt];
        elseif offset < 0
            currentParametersUpdated.k_on_mt     = currentParameters.k_on_mt(1:end+offset);
            currentParametersUpdated.k_step_mt   = currentParameters.k_step_mt((-offset+1):end);
            currentParametersUpdated.k_detach_mt = currentParameters.k_detach_mt((-offset+1):end);
        end
    end
    
    for i = 3:size(resultFolderCorrected, 1)
        currentParametersCorrected = updateParameters(currentParameters, resultFolderCorrected.offset(i) + correctionOffsetBinSizeLimit);
        simulateAndComputeSingleMeasurementModel(resultFolderCorrected.absFolderName{i}, parameterIndex, currentParametersCorrected, tspan);
    end
 end