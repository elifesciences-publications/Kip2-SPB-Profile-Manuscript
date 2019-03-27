function [fitDataXcoords, modelAtData, fitDataMean, fitDataSEM, currentBinCoords, currentBinProfile, shiftedLocations, meanGreenData, semGreenData] = computeModelAtData(currentBinCoords, currentBinProfile, laserFitParams, c )
    if all(isnan(currentBinCoords))
        error('Coords are NaN!');
    end
    
    currentBinProfile(~isfinite(currentBinProfile)) = 0;

    binCoordsAsc = currentBinCoords(2:end) - currentBinCoords(1:end-1);
    binCoordsOK = (isfinite(binCoordsAsc) & (binCoordsAsc > 0));

    while ~all(binCoordsOK)
        binCoordsOKconstrained = [false binCoordsOK];
        currentBinCoords = currentBinCoords(binCoordsOKconstrained);
        currentBinProfile = currentBinProfile(binCoordsOKconstrained);

        binCoordsAsc = currentBinCoords(2:end) - currentBinCoords(1:end-1);
        binCoordsOK = (isfinite(binCoordsAsc) & (binCoordsAsc > 0));
    end

    shiftedLocations        = c.lengthVectorCenterShifted;
    SPBlocation             = (c.meanMeanOffset-0.5)*(4/30);
    shiftedLocations        = -(shiftedLocations - SPBlocation);

    greenIntensitiesShifted = c.greenIntensitiesCenterShifted;
    bgIntensities = c.backgroundIntensities';

    if ~isempty(laserFitParams)
        greenIntensitiesShifted = laserFitParams(1) .* greenIntensitiesShifted + laserFitParams(2);
        bgIntensities = laserFitParams(1) .* bgIntensities + laserFitParams(2);
    end
    
    greenData = greenIntensitiesShifted - mean(bgIntensities);

    meanGreenData = mean(greenData, 2, 'omitnan');
    stdGreenData  = std(greenData, 0, 2, 'omitnan');
    nMTs = sum(~isnan(greenData),2);

    semGreenData = stdGreenData./sqrt(nMTs);

    [zero, zeroLoc] = min(abs(shiftedLocations));
    if zero > 1e-3
        error('SPB location failure');
    end

    fitDataIndices = (shiftedLocations >= -0.07 & shiftedLocations < (shiftedLocations(zeroLoc - 2*c.meanMeanOffset) + 0.07)); % + 0.07
    fitDataXcoords = shiftedLocations(fitDataIndices);

    fitDataMean = meanGreenData(fitDataIndices);
    fitDataSEM  = semGreenData(fitDataIndices);   

    lastAscendingIndex = find(currentBinCoords(2:end) <= currentBinCoords(1:end-1), 1, 'first');
    if ~isempty(lastAscendingIndex)
        warning('Parameter set with odd coordinates: %i', currentParameterSetIndex);
        currentBinCoords = currentBinCoords(1:lastAscendingIndex);
        currentBinProfile = currentBinProfile(1:lastAscendingIndex);
    end

    modelAtData = interp1qr(currentBinCoords', currentBinProfile', fitDataXcoords);
end

