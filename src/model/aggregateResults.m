function aggregateResults(i, resultFolderCorrected, absResultFolder, binsToCompute, exportPlots, debug, p)
    
    nComputedLengths = size(resultFolderCorrected, 1);
    nBins = size(binsToCompute, 1);
    %binProfiles = cell(nBins, 1);
    binIndex = 0;
    
    longestFile = load([resultFolderCorrected.absVmFolderName{nComputedLengths} int2str(i) '.mat']);
    %pxSize = mean(longestFile.result.model.meanCoords(2:end) - longestFile.result.model.meanCoords(1:end-1));
    binSize = longestFile.result.vm.params.binSize;
    
    maxLengthOverall = 0;
    maxIndOverall = 0;
    for j = 1:nComputedLengths
        currentFile = load([resultFolderCorrected.absVmFolderName{j} int2str(i) '.mat']);
        if currentFile.result.vm.binnedProfileMaxLengthOverall > maxLengthOverall
            maxLengthOverall = currentFile.result.vm.binnedProfileMaxLengthOverall;
        end
        currentIndOverall = 2*currentFile.result.vm.binnedProfileMaxLengthOverall + currentFile.result.vm.binnedPeakOverhangOverall - 1;
        if  currentIndOverall > maxIndOverall
            maxIndOverall = currentIndOverall;
        end       
    end
    
    maxIndAll = maxIndOverall; %- 1 

    nTimeSteps = size(longestFile.result.simResult.mtState, 1);
    
    result = struct;
    resultOut = struct;
    result.vmBin = cell(nBins, 1);
    resultOut.vmBin = cell(nBins, 1);
    
    


    for bin = binsToCompute'
        binIndex = binIndex + 1;
        nProfilesInBin = 0;
        
        totalNprofiles = nTimeSteps * longestFile.result.vm.params.binSize * nComputedLengths;
        result.vmBin{binIndex}.allAlignedBinnedProfiles = nan(totalNprofiles, maxIndAll);
        result.vmBin{binIndex}.allAlignedBinnedProfileStarts = nan(1, totalNprofiles);
        result.vmBin{binIndex}.allAlignedBinnedProfileCoords = nan(totalNprofiles, maxIndAll);
        result.vmBin{binIndex}.meanBoundKip2s = nan(1, nComputedLengths);
        result.vmBin{binIndex}.mtLengths = nan(1, totalNprofiles);
        
        for j = 1:nComputedLengths
            currentFile = load([resultFolderCorrected.absVmFolderName{j} int2str(i) '.mat']);   
            result.vmBin{binIndex}.meanBoundKip2s(j) = mean(sum(currentFile.result.simResult.mtState, 2));
            for fileJ = 1:binSize
                for fileI = 1:nTimeSteps
                    peakIndex = currentFile.result.vm.binnedPeakLocations(fileI, fileJ);
                    if isnan(peakIndex)
                        continue;
                    end
                    peakLocation = currentFile.result.vm.xCoordsBinned((peakIndex - 1)*binSize + fileJ);
                    
                    if length(bin) ~= 2

                        error(['Poop #' int2str(i) ': '  mat2str(binsToCompute)]); 
                    end
                    
                    if (peakLocation < bin(1) || peakLocation >  bin(2))
                        continue;
                    end
                    nProfilesInBin = nProfilesInBin + 1;
                    % Peak location is within current bin
                    startIndMinus1 = maxLengthOverall - peakIndex;

                    % Index needs to consider MT length
                    profileIndex = (j-1)*binSize*nTimeSteps + (fileJ-1)*nTimeSteps + fileI;
                    
                    result.vmBin{binIndex}.allAlignedBinnedProfiles(profileIndex, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 1)) = currentFile.result.vm.binnedProfiles{fileJ}(fileI, :);
%                     if startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 1 > maxIndAll
%                         warning('moo');
%                         startIndMinus1 + 1
%                         startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 1
%                     end
                    result.vmBin{binIndex}.allAlignedBinnedProfiles(profileIndex, (startIndMinus1 + 2):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 2)) = 0.5 * (...
                        result.vmBin{binIndex}.allAlignedBinnedProfiles(profileIndex, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 3)) + ...
                        result.vmBin{binIndex}.allAlignedBinnedProfiles(profileIndex, (startIndMinus1 + 3):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 1)) ...
                    );                    
                    result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 1)) = currentFile.result.vm.binnedXCoords{fileJ};

                    result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, (startIndMinus1 + 2):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 2)) = 0.5 * (...
                        result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 3)) + ...
                        result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, (startIndMinus1 + 3):2:(startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :)) - 1)) ...
                    );               
                    for k = startIndMinus1:-1:1
                        result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, k) = result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, k+1) - currentFile.result.vm.params.binSize*0.008 * 0.5;
                    end

                    for k = (startIndMinus1 + 2*length(currentFile.result.vm.binnedProfiles{fileJ}(fileI, :))):maxIndAll
                        result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex,k) = result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, k-1) + currentFile.result.vm.params.binSize*0.008 * 0.5;
                    end           

                    %result.vmBin{binIndex}.allAlignedBinnedProfiles(profileIndex, :)
                    %result.vmBin{binIndex}.allAlignedBinnedProfileCoords(profileIndex, :)
                    result.vmBin{binIndex}.allAlignedBinnedProfileStarts(profileIndex) = startIndMinus1 + 1;
                    result.vmBin{binIndex}.mtLengths(profileIndex) = currentFile.result.simResult.parameters.maxLength;
                    
                    % result.vmBin{binIndex}.mean
                end
            end
        end
        

        
        resultOut.vmBin{binIndex}.meanAllAlignedBinnedProfileCoords = mean(result.vmBin{binIndex}.allAlignedBinnedProfileCoords, 1, 'omitnan');
        resultOut.vmBin{binIndex}.meanAllAlignedBinnedProfile       = mean(result.vmBin{binIndex}.allAlignedBinnedProfiles, 1, 'omitnan');
        resultOut.vmBin{binIndex}.meanAllAlignedBinnedProfileN      = sum(~isnan(result.vmBin{binIndex}.allAlignedBinnedProfiles));
        resultOut.vmBin{binIndex}.nProfilesInBin = nProfilesInBin;
        resultOut.vmBin{binIndex}.meanBoundKip2s = result.vmBin{binIndex}.meanBoundKip2s;
    end
    
    if (exportPlots || debug)
        fig = figure();
        hold on;
        colors = get(gca,'colororder');
        binIndex = 0;
        for bin = binsToCompute'
            binIndex = binIndex + 1;
            nonNanIndices = find(~isnan(result.vmBin{binIndex}.allAlignedBinnedProfileStarts));
            rawCurvesToPlot = nonNanIndices(randsample(length(nonNanIndices), 100)); 
            plot(result.vmBin{binIndex}.allAlignedBinnedProfileCoords(rawCurvesToPlot, :)', result.vmBin{binIndex}.allAlignedBinnedProfiles(rawCurvesToPlot, :)', 'Color', [colors(binIndex, :) 0.1])
        end
        binIndex = 0;
        for bin = binsToCompute'
            binIndex = binIndex + 1;        
            h = plot(resultOut.vmBin{binIndex}.meanAllAlignedBinnedProfileCoords, resultOut.vmBin{binIndex}.meanAllAlignedBinnedProfile, 'Color', colors(binIndex, :), 'LineWidth', 1.5);
        end
        plot([1 1] .* mean(result.vmBin{binIndex}.mtLengths, 'omitnan') .* 0.008, [0 1], '--k', 'LineWidth', 1.5);
        
        legend([h], {sprintf('nKip2Free = %g nM\nk_{on} = %g\nk_{in} = %g\nk_{out} = %g', p(1:4))});
        legend boxoff
        ylim([0 1]);
        xlim([-0.5 4]);
        xlabel('Distance from SPB ({\mu}m)');
        ylabel('Occupancy');
        applyPaperFormatting();
        if exportPlots
            %print([resultFolderCorrectedMeanVmImages filesep int2str(i) '.pdf'], '-dpdf');
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
    
    result = resultOut;
    result.binsToCompute = binsToCompute;
    parsave([absResultFolder.correctedMeanVm int2str(i) '.mat'], result);
end