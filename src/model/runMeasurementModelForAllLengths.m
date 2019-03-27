function runMeasurementModelForAllLengths(resultFolderCorrected, parameterIndex, currentParameters, tspan)
    % © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)
    function currentParametersUpdated = updateParameters(currentParameters, mtLength)
        currentParametersUpdated = currentParameters;
        currentParametersUpdated.maxLength = mtLength;
        
        mtLengthCurrent = length(currentParameters.k_step_mt);
        if mtLength >= mtLengthCurrent
            currentParametersUpdated.k_on_mt     = [currentParameters.k_on_mt currentParameters.k_on_mt(end) * ones(1, mtLength - mtLengthCurrent)];
            currentParametersUpdated.k_step_mt   = [currentParameters.k_step_mt(1) * ones(1, mtLength - mtLengthCurrent) currentParameters.k_step_mt];
            currentParametersUpdated.k_detach_mt = [currentParameters.k_detach_mt(1) * ones(1, mtLength - mtLengthCurrent) currentParameters.k_detach_mt];
        elseif mtLength < length(currentParameters.k_step_mt)
            currentParametersUpdated.k_on_mt     = currentParameters.k_on_mt(1:mtLength);
            currentParametersUpdated.k_step_mt   = currentParameters.k_step_mt((end-mtLength+1):end);
            currentParametersUpdated.k_detach_mt = currentParameters.k_detach_mt((end-mtLength+1):end);
        end
    end
    
    szResultFolder = size(resultFolderCorrected, 1);
    for i = 1:szResultFolder
        currentParametersCorrected = updateParameters(currentParameters, currentParameters.lenghtsToCompute(i));
        runMeasurementModel(resultFolderCorrected.absFolderName{i}, parameterIndex, currentParametersCorrected, tspan);
    end
 end