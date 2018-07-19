function conditionResultOut = createTwoColorSubset(conditionResult, indices, newName)
    conditionResultOut = struct;
    conditionResultOut.totalLengths                  = conditionResult.totalLengths(indices);
    conditionResultOut.totalLengthsUM                = conditionResult.totalLengthsUM(indices);
    conditionResultOut.lengthsUM                     = conditionResult.lengthsUM(indices);
    conditionResultOut.lengths                       = conditionResult.lengths(indices);
    conditionResultOut.condition                     = newName;
    
    
    conditionResultOut.lengthVector                  = conditionResult.lengthVector;
    conditionResultOut.greenIntensities              = conditionResult.greenIntensities(:, indices);
    conditionResultOut.redIntensities                = conditionResult.redIntensities(:, indices);
    conditionResultOut.offset                        = conditionResult.offset;
    conditionResultOut.offsetSPB                     = conditionResult.offsetSPB;
    conditionResultOut.greenIntensitiesShifted       = conditionResult.greenIntensitiesShifted(:, indices);
    conditionResultOut.redIntensitiesShifted         = conditionResult.redIntensitiesShifted(:, indices);
    conditionResultOut.lengthVectorShifted           = conditionResult.lengthVectorShifted;
    conditionResultOut.greenIntensitiesSPBShifted    = conditionResult.greenIntensitiesSPBShifted(:, indices);
    conditionResultOut.redIntensitiesSPBShifted      = conditionResult.redIntensitiesSPBShifted(:, indices);
    conditionResultOut.lengthVectorSPBShifted        = conditionResult.lengthVectorSPBShifted;
    
    conditionResultOut.lengthVectorCenterShifted     = conditionResult.lengthVectorCenterShifted;
    conditionResultOut.greenIntensitiesCenterShifted = conditionResult.greenIntensitiesCenterShifted;
    conditionResultOut.redIntensitiesCenterShifted   = conditionResult.redIntensitiesCenterShifted;
    
    conditionResultOut.greenPeaks                    = conditionResult.greenPeaks(indices);
    conditionResultOut.redPeaks                      = conditionResult.redPeaks(indices);
    
    conditionResultOut.cellNames                     = conditionResult.cellNames(indices);
end

