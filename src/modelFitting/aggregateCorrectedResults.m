function aggregateCorrectedResults(i, resultFolderCorrected, absResultFolder, zeroOffsetIndex, exportPlots, debug, parameterCombinations)
    
    zeroOffsetFile = load([resultFolderCorrected.absVmFolderName{zeroOffsetIndex} int2str(i) '.mat']);
    
    if zeroOffsetFile.result.model.nProfiles > 0
        pxSize = mean(zeroOffsetFile.result.model.meanCoords(2:end) - zeroOffsetFile.result.model.meanCoords(1:end-1));

        zeroOffsetPeakLocation = zeroOffsetFile.result.model.meanCoords(end) - zeroOffsetFile.result.model.meanPlusEndLengthOffset;

        leftBound  = zeroOffsetFile.result.model.meanCoords(1);
        rightBound = zeroOffsetFile.result.model.meanCoords(end);
        for j = 3:size(resultFolderCorrected, 1)
            currentFile = load([resultFolderCorrected.absVmFolderName{j} int2str(i) '.mat']);

            if ~isempty(currentFile.result.model.meanCoords)
                currentOffsetPeakLocation = currentFile.result.model.meanCoords(end) - currentFile.result.model.meanPlusEndLengthOffset;

                alignmentOffset = (zeroOffsetPeakLocation - currentOffsetPeakLocation) * 0.5;

                currentLeftBound = zeroOffsetFile.result.model.meanCoords(1) + alignmentOffset;
                currentRightBound = zeroOffsetFile.result.model.meanCoords(end) + alignmentOffset;

                leftBound  = min([leftBound,  currentLeftBound] , [], 'omitnan');
                rightBound = max([rightBound, currentRightBound], [], 'omitnan');
            end
        end

        totalLength = rightBound - leftBound;
        interpLength = round(totalLength/pxSize * 2); % twice the resolution of the individual results
        xCoords = linspace(leftBound, rightBound, interpLength + 1); 

        intensities = nan(size(resultFolderCorrected, 1) - 2, length(xCoords));
        Ndata =  nan(size(resultFolderCorrected, 1) - 2, length(xCoords));
        weights = nan(size(resultFolderCorrected, 1) - 2, 1);
        correctedLength = NaN;
        for j = 3:size(resultFolderCorrected, 1)
            currentFile = load([resultFolderCorrected.absVmFolderName{j} int2str(i) '.mat']);

            weights(j-2) = currentFile.result.model.nProfiles;
            if currentFile.result.model.nProfiles > 0
                currentOffsetPeakLocation = currentFile.result.model.meanCoords(end) - currentFile.result.model.meanPlusEndLengthOffset;

                alignmentOffset = (zeroOffsetPeakLocation - currentOffsetPeakLocation) * 0.5;
                intensities(j-2, :) = interp1(currentFile.result.model.meanCoords + alignmentOffset, currentFile.result.model.mean, xCoords);
                Ndata(j-2, :) = interp1(currentFile.result.model.meanCoords + alignmentOffset, currentFile.result.model.meanN, xCoords);
            end

            if resultFolderCorrected.offset(j) == 0
                correctedLength = size(currentFile.result.simResult.mtState, 2);
            end
        end
        result = struct();
        result.xCoords = xCoords;
        result.intensities = intensities;
        result.Ndata = Ndata;
        result.weights = weights;

        Ndata(isnan(Ndata)) = 0;
        NdataSum = sum(Ndata, 1);
        NdataSumInv = 1 ./ NdataSum;

        tempIntensities = intensities;
        tempIntensities(isnan(tempIntensities)) = 0;

        result.weightedIntensities = sum((tempIntensities .* Ndata), 1) .* NdataSumInv;
    else
        result = struct();
        result.xCoords = [];
        result.intensities = [];
        result.Ndata = [];
        result.weights = [];

        result.weightedIntensities = [];        
    end
    

    

    
    if (exportPlots || debug) && zeroOffsetFile.result.model.nProfiles > 0
        fig = figure();
        hold on;
        h2 = plot(xCoords, intensities);
        tempIntensities = intensities';
        tempIntensities(isnan(tempIntensities)) = 0;
        
        totalWeightInv = 1/sum(weights);
        plot(xCoords, mean(intensities, 'omitnan'), '-r', 'LineWidth', 2);
        plot(xCoords, (tempIntensities * weights .* totalWeightInv)' , '-k', 'LineWidth', 2);
        
        h = plot(xCoords, result.weightedIntensities  , '-b', 'LineWidth', 2);
        
        p = parameterCombinations(:, i);
        legend([h; h2], [{sprintf('nKip2Free = %g nM\nk_{on} = %g\nk_{in} = %g\nk_{out} = %g', p(2:5))}, arrayfun(@int2str, resultFolderCorrected.offset(3:end)', 'UniformOutput', false)]);
        legend boxoff
        title(sprintf('Raw length = %i, corrected length = %i', p(1), correctedLength));
        ylim([0 1]);
        xlim([-1 3]);
        xlabel('Distance from SPB ({\mu}m)');
        ylabel('Occupancy');
        if exportPlots
            done = 0;
            while done >= 0 && done < 10
                try
                    print([absResultFolder.correctedMeanVmImages int2str(i) '.png'], '-dpng');
                    done = -1;
                catch ME
                    fprintf(['Error exporting figure: ' ME.message '\n']);
                    pause(1);
                    done = done + 1;
                end
            end
            if done >= 10
                error('Failure');
            end
        end
        if ~debug
            close(fig);
        end
    end
    
    parsave([absResultFolder.correctedMeanVm int2str(i) '.mat'], result);
end