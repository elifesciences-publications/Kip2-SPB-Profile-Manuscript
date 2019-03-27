function [fitresult, gof] = laserPowerModel(cFrom, cRef)
    quantilesForCalibration = [0.005:0.002:0.995];
    nMicrotubulesRef = size(cRef.greenIntensities, 2);
    refIntensities = [];

    refPlusEndLoc = cRef.offset;

    for j = 1:nMicrotubulesRef
        if cRef.lengths(j) ~= 0
            refSPBloc = refPlusEndLoc + cRef.lengths(j);
            currentMtIntensities = cRef.greenIntensitiesShifted(refPlusEndLoc:refSPBloc, j);
            refIntensities = [refIntensities; currentMtIntensities];
        end
    end

    assert(all(isfinite(refIntensities)));

    nMicrotubulesFrom = size(cFrom.greenIntensities, 2);
    fromIntensities = [];

    fromPlusEndLoc = cFrom.offset;

    for j = 1:nMicrotubulesFrom
        if cFrom.lengths(j) ~= 0
            fromSPBloc = fromPlusEndLoc + cFrom.lengths(j);
            currentMtIntensities = cFrom.greenIntensitiesShifted(fromPlusEndLoc:fromSPBloc, j);
            fromIntensities = [fromIntensities; currentMtIntensities];
        end
    end            

    assert(all(isfinite(fromIntensities))); 

    allQuantiles = quantile(fromIntensities, quantilesForCalibration);
    allQuantilesRef = quantile(refIntensities, quantilesForCalibration);

    [fitresult, gof] = laserModelFit(allQuantiles, allQuantilesRef, false);
end