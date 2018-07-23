function result = computeMeasurementModel(simResult, outFileName, binBounds)
    %% 1. Perform virtual microscopy
    
    mode = 'both'; % others: SPB, both
    useBinBounds = exist('binBounds', 'var');  
    vm = struct;
    
    p = struct;
    
    lambda_em = 509;
    NA = 1.46;
    sigma_im = 0.21 * lambda_em / NA;
    
    p.padding = 40; % padding past MT
    p.lambda_em = lambda_em; % fluorophore emission wavelength in nanometers
    p.NA = NA; % numerical aperture
    p.sigma_im = sigma_im; % gaussian sigma in nanometers
    p.sigma_sites = sigma_im / 8; % gaussian sigma in 8 nm kinesin sites
    p.binSize = 17; % pixel size in 8 nm kinesin sites.
    p.debug = true;
    vm.params = p;

    %sites_per_pixel = 17; % 8 nm * 17 = 136 nm (microscope: 133.1 nm)


    mtState = simResult.mtState;

    nTimeSteps = size(mtState, 1);
    curLength = size(mtState, 2);

    result = struct;
    result.simResult = simResult;
    result.vm = vm;
    
    % 1.1 Convolve with PSF
    result.vm.paddedState = [zeros(nTimeSteps, vm.params.padding) mtState zeros(nTimeSteps, vm.params.padding)];
    result.vm.convMat = imgaussfilt(result.vm.paddedState, [1e-10 vm.params.sigma_sites]); % use very low spread in time
    result.vm.xCoordsConv = ((1:size(result.vm.convMat,2))-1-result.vm.params.padding) .* 0.008;  
    
    % 1.2 Compute all possible binnings
    convMatBinPadded = [zeros(nTimeSteps, (vm.params.binSize - 1) / 2) result.vm.convMat zeros(nTimeSteps, (vm.params.binSize - 1) / 2)];
    result.vm.convMatBinned = nan(nTimeSteps, size(convMatBinPadded, 2));
    maxSize = size(convMatBinPadded, 2);
    for i = 1:maxSize
        minIndex = max(1, i - (vm.params.binSize - 1) / 2);
        maxIndex = min(maxSize, i + (vm.params.binSize - 1) / 2 );
        result.vm.convMatBinned(:, i) = mean(convMatBinPadded(:, minIndex:maxIndex), 2);
    end
    result.vm.xCoordsBinned = ((1:maxSize) - 1 - vm.params.padding - (vm.params.binSize - 1)/2)*0.008;
    
    % 1.3 Sample all possible binnings
    result.vm.binnedProfiles = {};
    result.vm.binnedProfileMaxLength = zeros(1, vm.params.binSize);
    result.vm.binnedXCoords = {};
    for i = 1:vm.params.binSize
        result.vm.binnedProfiles{i} = result.vm.convMatBinned(:, i:vm.params.binSize:end);
        result.vm.binnedProfileMaxLength(i) = max(size(result.vm.binnedProfiles{i}, 2), result.vm.binnedProfileMaxLength(i));
        result.vm.binnedXCoords{i} = result.vm.xCoordsBinned(i:vm.params.binSize:end);
    end
    
    result.vm.binnedProfileMaxLengthOverall = max(result.vm.binnedProfileMaxLength);
    
    % 1.4 Do peak detection
    result.vm.binnedPeakLocations = nan(nTimeSteps, vm.params.binSize);
    result.vm.binnedPeakOverhang = zeros(1, vm.params.binSize);

    for j = 1:vm.params.binSize
        result.vm.hasPeak{j} = zeros(1, nTimeSteps);
        result.vm.isInBin{j} = zeros(1, nTimeSteps);
        for i = 1:nTimeSteps
            greenResult = struct;
            [~, greenResult.locs] = findpeaks_stripped(result.vm.binnedProfiles{j}(i, :)); % optimize for speed
            % [greenResult.pks, greenResult.locs, greenResult.w, greenResult.p] = findpeaks(result.vm.binnedProfiles{j}(i, :));
            
            % Plus end peak is the last one, if it exists
            result.vm.hasPeak{j}(i) = ~isempty(greenResult.locs);
            
            if ~result.vm.hasPeak{j}(i)
                result.vm.isInBin{j}(i) = false;
            elseif useBinBounds 
                peakIndex = greenResult.locs(end);
                peakLocation = result.vm.xCoordsBinned((peakIndex - 1)*vm.params.binSize + j);
                %peakLocation = 17*round((peakLocation/0.008)/17)*0.008; % round to nearest bin
                result.vm.isInBin{j}(i) = (binBounds(1) <= peakLocation && binBounds(2) >= peakLocation);
            else
                result.vm.isInBin{j}(i) = true;
            end
           
            
            if result.vm.isInBin{j}(i)
                result.vm.binnedPeakLocations(i, j) = greenResult.locs(end);
            end
                
            % SPB is assumed to be at 0, so no peak detection here
            
        end
        
        result.vm.binnedPeakOverhang(j) = result.vm.binnedProfileMaxLength(j) - min(result.vm.binnedPeakLocations(:, j), [], 'omitnan');
    end
    result.vm.binnedPeakOverhangOverall = max(result.vm.binnedPeakOverhang, [], 'omitnan');

    
    % 1.5 Align profiles
    result.vm.alignedBinnedProfiles = cell(1, vm.params.binSize);
    result.vm.alignedBinnedProfileStarts = nan(nTimeSteps, vm.params.binSize);
    result.vm.alignedBinnedProfileCoords = cell(1, vm.params.binSize);
    
    switch mode
        case 'plusEnd'
            maxIndAll = result.vm.binnedProfileMaxLengthOverall + result.vm.binnedPeakOverhangOverall;
            if isnan(maxIndAll)
                maxIndAll = 0;
            end
            result.vm.allAlignedBinnedProfiles = nan(nTimeSteps * vm.params.binSize, maxIndAll);
            result.vm.allAlignedBinnedProfileStarts = nan(1, nTimeSteps * vm.params.binSize);
            result.vm.allAlignedBinnedProfileCoords = nan(nTimeSteps * vm.params.binSize, maxIndAll);

            for j = 1:vm.params.binSize
                maxInd = result.vm.binnedProfileMaxLength(j) + result.vm.binnedPeakOverhang(j); 

                if isnan(maxInd)
                    maxInd = 0;
                end
                result.vm.alignedBinnedProfiles{j} = nan(nTimeSteps, maxInd);
                result.vm.alignedBinnedProfileCoords{j} = nan(nTimeSteps, maxInd);

                for i = 1:nTimeSteps
                    if ~result.vm.isInBin{j}(i) 
                        continue;
                    end

                    % Plus end peak is the last one
                    startIndMinus1 = result.vm.binnedProfileMaxLength(j) - result.vm.binnedPeakLocations(i,j);

                    result.vm.alignedBinnedProfiles{j}(i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedProfiles{j}(i, :);
                    result.vm.alignedBinnedProfileCoords{j}(i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedXCoords{j};
                    for k = startIndMinus1:-1:1
                        result.vm.alignedBinnedProfileCoords{j}(i,k) = result.vm.alignedBinnedProfileCoords{j}(i,k+1) - vm.params.binSize*0.008;
                    end

                    for k = (startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)) + 1):maxInd
                        result.vm.alignedBinnedProfileCoords{j}(i,k) = result.vm.alignedBinnedProfileCoords{j}(i,k-1) + vm.params.binSize*0.008;
                    end
                    result.vm.alignedBinnedProfileStarts(i, j) = startIndMinus1 + 1;


                    startIndMinus1 = result.vm.binnedProfileMaxLengthOverall - result.vm.binnedPeakLocations(i,j);
                    result.vm.allAlignedBinnedProfiles((j-1)*nTimeSteps + i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedProfiles{j}(i, :);
                    result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedXCoords{j};

                    for k = startIndMinus1:-1:1
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i,k) = result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, k+1) - vm.params.binSize*0.008;
                    end

                    for k = (startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)) + 1):maxIndAll
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i,k) = result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, k-1) + vm.params.binSize*0.008;
                    end           


                    result.vm.allAlignedBinnedProfileStarts((j-1)*nTimeSteps + i) = startIndMinus1 + 1;
                end
            end



        case 'SPB'
            maxIndAll = result.vm.binnedProfileMaxLengthOverall;
            if isnan(maxIndAll)
                maxIndAll = 0;
            end
            result.vm.allAlignedBinnedProfiles = nan(nTimeSteps * vm.params.binSize, maxIndAll);
            result.vm.allAlignedBinnedProfileStarts = nan(1, nTimeSteps * vm.params.binSize);
            result.vm.allAlignedBinnedProfileCoords = nan(nTimeSteps * vm.params.binSize, maxIndAll);
            
            for j = 1:vm.params.binSize
                maxInd = result.vm.binnedProfileMaxLength(j); 

                if isnan(maxInd)
                    maxInd = 0;
                end
                result.vm.alignedBinnedProfiles{j} = nan(nTimeSteps, maxInd);
                result.vm.alignedBinnedProfileCoords{j} = nan(nTimeSteps, maxInd);

                for i = 1:nTimeSteps
                    if ~result.vm.isInBin{j}(i) 
                        continue;
                    end

                    % Minus end aligned
                    startIndMinus1 = 0;

                    result.vm.alignedBinnedProfiles{j}(i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedProfiles{j}(i, :);
                    result.vm.alignedBinnedProfileCoords{j}(i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedXCoords{j};
                    for k = startIndMinus1:-1:1
                        result.vm.alignedBinnedProfileCoords{j}(i,k) = result.vm.alignedBinnedProfileCoords{j}(i,k+1) - vm.params.binSize*0.008;
                    end

                    for k = (startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)) + 1):maxInd
                        result.vm.alignedBinnedProfileCoords{j}(i,k) = result.vm.alignedBinnedProfileCoords{j}(i,k-1) + vm.params.binSize*0.008;
                    end
                    result.vm.alignedBinnedProfileStarts(i, j) = startIndMinus1 + 1;


                    startIndMinus1 = 0;
                    result.vm.allAlignedBinnedProfiles((j-1)*nTimeSteps + i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedProfiles{j}(i, :);
                    result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, (startIndMinus1 + 1):(startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)))) = result.vm.binnedXCoords{j};

                    for k = startIndMinus1:-1:1
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i,k) = result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, k+1) - vm.params.binSize*0.008;
                    end

                    for k = (startIndMinus1 + length(result.vm.binnedProfiles{j}(i, :)) + 1):maxIndAll
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i,k) = result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, k-1) + vm.params.binSize*0.008;
                    end           


                    result.vm.allAlignedBinnedProfileStarts((j-1)*nTimeSteps + i) = startIndMinus1 + 1;

                    % SPB is assumed to be at 0, so no peak detection here

                end
            end
            
            
        case 'both'
            maxIndAll = 2*result.vm.binnedProfileMaxLengthOverall - 1 + result.vm.binnedPeakOverhangOverall;
            if isnan(maxIndAll)
                maxIndAll = 0;
            end
            result.vm.allAlignedBinnedProfiles = nan(nTimeSteps * vm.params.binSize, maxIndAll);
            result.vm.allAlignedBinnedProfileStarts = nan(1, nTimeSteps * vm.params.binSize);
            result.vm.allAlignedBinnedProfileCoords = nan(nTimeSteps * vm.params.binSize, maxIndAll);

            for j = 1:vm.params.binSize
                maxInd = 2*result.vm.binnedProfileMaxLength(j) - 1 + result.vm.binnedPeakOverhang(j); 

                if isnan(maxInd)
                    maxInd = 0;
                end
                result.vm.alignedBinnedProfiles{j} = nan(nTimeSteps, maxInd);
                result.vm.alignedBinnedProfileCoords{j} = nan(nTimeSteps, maxInd);

                for i = 1:nTimeSteps
                    if ~result.vm.isInBin{j}(i) 
                        continue;
                    end

                    % Plus end peak is the last one
                    startIndMinus1 = result.vm.binnedProfileMaxLength(j) - result.vm.binnedPeakLocations(i,j);

                    result.vm.alignedBinnedProfiles{j}(i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) = result.vm.binnedProfiles{j}(i, :);
                    result.vm.alignedBinnedProfiles{j}(i, (startIndMinus1 + 2):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 2)) = 0.5 * (...
                        result.vm.alignedBinnedProfiles{j}(i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 3)) + ...
                        result.vm.alignedBinnedProfiles{j}(i, (startIndMinus1 + 3):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) ...
                    );
                    result.vm.alignedBinnedProfileCoords{j}(i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) = result.vm.binnedXCoords{j};
                    
                    result.vm.alignedBinnedProfileCoords{j}(i, (startIndMinus1 + 2):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 2)) = 0.5 * (...
                        result.vm.alignedBinnedProfileCoords{j}(i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 3)) + ...
                        result.vm.alignedBinnedProfileCoords{j}(i, (startIndMinus1 + 3):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) ...
                    );
                    
                    
                    for k = startIndMinus1:-1:1
                        result.vm.alignedBinnedProfileCoords{j}(i,k) = result.vm.alignedBinnedProfileCoords{j}(i,k+1) - vm.params.binSize*0.008 * 0.5;
                    end

                    for k = (startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :))):maxInd
                        result.vm.alignedBinnedProfileCoords{j}(i,k) = result.vm.alignedBinnedProfileCoords{j}(i,k-1) + vm.params.binSize*0.008 * 0.5;
                    end
                    result.vm.alignedBinnedProfileStarts(i, j) = startIndMinus1 + 1;


                    startIndMinus1 = result.vm.binnedProfileMaxLengthOverall - result.vm.binnedPeakLocations(i,j);
                    result.vm.allAlignedBinnedProfiles((j-1)*nTimeSteps + i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) = result.vm.binnedProfiles{j}(i, :);
                    
                    result.vm.allAlignedBinnedProfiles((j-1)*nTimeSteps + i, (startIndMinus1 + 2):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 2)) = 0.5 * (...
                        result.vm.allAlignedBinnedProfiles((j-1)*nTimeSteps + i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 3)) + ...
                        result.vm.allAlignedBinnedProfiles((j-1)*nTimeSteps + i, (startIndMinus1 + 3):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) ...
                    );                    
                    
                    result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) = result.vm.binnedXCoords{j};

                    result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, (startIndMinus1 + 2):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 2)) = 0.5 * (...
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, (startIndMinus1 + 1):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 3)) + ...
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, (startIndMinus1 + 3):2:(startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :)) - 1)) ...
                    );               
                    
                    for k = startIndMinus1:-1:1
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i,k) = result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, k+1) - vm.params.binSize*0.008 * 0.5;
                    end

                    for k = (startIndMinus1 + 2*length(result.vm.binnedProfiles{j}(i, :))):maxIndAll
                        result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i,k) = result.vm.allAlignedBinnedProfileCoords((j-1)*nTimeSteps + i, k-1) + vm.params.binSize*0.008 * 0.5;
                    end           


                    result.vm.allAlignedBinnedProfileStarts((j-1)*nTimeSteps + i) = startIndMinus1 + 1;
                end
            end            
            
        case default
            error(['Unknown mode: ' mode]);
    end
    result.vm.meanAlignedBinnedProfileCoords = cell(1, vm.params.binSize);
    for j = 1:vm.params.binSize          
        result.vm.meanAlignedBinnedProfileCoords{j} = mean(result.vm.alignedBinnedProfileCoords{j},1, 'omitnan');
    end

    result.vm.meanAllAlignedBinnedProfileCoords = mean(result.vm.allAlignedBinnedProfileCoords,1, 'omitnan');
    %% 2. Average profiles
    result.model.mean = mean(result.vm.allAlignedBinnedProfiles(:,:), 1, 'omitnan');
    result.model.meanN = sum(~isnan(result.vm.allAlignedBinnedProfiles(:,:)), 1);
    result.model.meanCoords = result.vm.meanAllAlignedBinnedProfileCoords;
    
    %% 3. Compute peak to true plus end offset and minus end offset
    result.model.plusEndLengthOffsets  = nan(nTimeSteps, vm.params.binSize);
    

    %result.model.plusEndLengthOffsetsAlt = cell(1, vm.params.binSize);
    for i = 1:vm.params.binSize
        % MT length - (length from start start to peak)
        %result.model.plusEndLengthOffsetsAlt{i}  = curLength*0.008 - ((result.vm.binnedPeakLocations{i} - 1)*vm.params.binSize - vm.params.padding - (vm.params.binSize - 1)/2 + i - 1)*0.008; % in sites
        hasPeak = isfinite(result.vm.binnedPeakLocations(:, i));
        result.model.plusEndLengthOffsets(hasPeak, i) = curLength*0.008 - result.vm.xCoordsBinned((result.vm.binnedPeakLocations(hasPeak, i) - 1)*vm.params.binSize + i);
        
    end 
    
    meanOmitNan = @(x) mean(x, 'omitnan');
    result.model.meanPlusEndLengthOffset = meanOmitNan(result.model.plusEndLengthOffsets(:));
    result.model.nProfiles = sum(isfinite(result.vm.allAlignedBinnedProfileStarts));
    
    result.model.meanMinusEndLengthOffset = nan;
    if ~isempty(result.model.meanCoords)
        result.model.meanMinusEndLengthOffset = -result.model.meanCoords(1);
    end
    
    parsave(outFileName, result);
end