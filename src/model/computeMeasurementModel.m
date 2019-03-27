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
            Yin = result.vm.binnedProfiles{j}(i, :);
            
            greenResult.locs = findpeaks_fast(Yin);
            
            % Plus end peak is the last one, if it exists
            result.vm.hasPeak{j}(i) = ~isempty(greenResult.locs);
            
            if ~result.vm.hasPeak{j}(i)
                result.vm.isInBin{j}(i) = false;
            elseif useBinBounds 
                peakIndex = greenResult.locs(end);
                peakLocation = result.vm.xCoordsBinned((peakIndex - 1)*vm.params.binSize + j);
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
    parsave(outFileName, result);
end