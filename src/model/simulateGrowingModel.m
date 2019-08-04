function result = simulateGrowingModel(parameters, tspan, currentFileName)
% simulateModel: Stochastic simulation of Kip2 motor 
%                model using the Direct Gillespie Algorithm
% 
% Arguments:
%   parameters: struct
%     .maxLength: double
%       length of filament
%     .k_on_mt: double vector
%       on rate for each site (nM^-1 s^-1)
%     .k_step_mt: double vector
%       step rate for each site (s^-1)
%     .k_detach_mt: double vector
%       off rate for each site (s^-1)
%     .nKip2Free: int
%       number of free motors at t = 0
%     .reportRunLengths: boolean
%       keep track of motor run lengths / run times
%       (turn off for performance reasons)
%   parameterIndex: double
%     index of parameter set (returned as part of the result struct)
%   tspan: double vector
%     vector of time points to export
%   currentFileName: char vector
%     file in which to save the result struct
%
%currentParameters.startLength      = lengthsToCompute;
%currentParameters.endLength        = lengthsToCompute(1);
%currentParameters.warmupTime       = warmupTime;
%currentParameters.growthRate       = mtGrowthRate;

%
%
% Output:
%   result: struct
%     .mtState: double matrix
%       matrix of length(tspan) x maxLength dimension containing
%       all requested samples (entry = 0: free site, 1: Kip2 bound)
%     .runLengthsKip2: double vector
%       motor run lengths (from binding to unbinding)
%     .runTimesKip2: double vector
%       motor run times (from binding to unbinding)
%     .steps: double
%       number of simulation steps taken
%     .runtime: double
%       runtime of algorithm in seconds
%     .tspan: double vector
%       vector of time points to export
%     .parameterIndex: double
%       index of parameter set specified when calling this function

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

    localStart = tic;

    %% Parameters  
    startLength       = parameters.startLength;
    endLength         = parameters.endLength;
    warmupTime        = parameters.warmupTime;
    
    k_on_mt           = parameters.k_on_mt;
    k_on_in_mt        = parameters.k_on_in_mt;
    k_step_mt         = parameters.k_step_mt;
    k_detach_mt       = parameters.k_detach_mt;
    k_detach_mt_end   = parameters.k_detach_mt_end;
    
    nKip2Free         = parameters.nKip2Free;
    
    growthRate        = parameters.growthRate;
    
    
    tend = tspan(end); % seconds
    mt = nan(1, endLength);
    mt(1:startLength) = 0;
    

    nReactions = 4;
    a = zeros(1, nReactions);
    
    step = 0;
    
    t = 0;
    tIndex = 1;
    mtState = nan(length(tspan), endLength);
    kip2FreeState = zeros(length(tspan), 1);
    
    movingKip2s = mt;
    tState = zeros(size(tspan));
    
    currentMtLength = startLength;
    while t < tend && currentMtLength <= endLength
        
        if nKip2Free < 0 || nKip2Free > parameters.nKip2Free
            error('Molecules not conserved!');
        end
        
        % nan: no microtubule
        % 0: Free
        % 1: Kip2
        
        while tIndex <= length(tspan) && t > tspan(tIndex)
            mtState(tIndex, :) = mt;
            tState(tIndex) = t;
            kip2FreeState(tIndex) = nKip2Free;
            tIndex = tIndex + 1;
        end

        %% Attachment               
        r_on_mt = k_on_mt * (60 / 140) * nKip2Free;
        r_on_in_mt = k_on_in_mt * (60 / 140) * nKip2Free;
        
        attachmentRates = double(~mt(1:currentMtLength)) * r_on_mt; % on the microtubule lattice
        attachmentRates(1) =  double(~mt(1)) * r_on_in_mt; % at the minus end
        
        attachmentCumSum = cumsum(attachmentRates);
        
        a(:, 1) = attachmentCumSum(end);
        
        %% Stepping
        stepRates = double(movingKip2s(1:currentMtLength)) * k_step_mt;
        stepRates(currentMtLength) = 0; % No stepping at the end of the microtubule
        stepCumSum = cumsum(stepRates);
        
        a(:, 2) = stepCumSum(end);
        
        %% Detachment   
        detachmentRates = double(mt(1:currentMtLength)) * k_detach_mt;
        detachmentRates(currentMtLength) = mt(currentMtLength) * k_detach_mt_end;
        detachmentCumSum = cumsum(detachmentRates);
        a(:, 3) = detachmentCumSum(end);    
        
        %% MT Growth Rate
        mtGrowthRateTotal = (t >= warmupTime) * growthRate;
        a(:, 4) = mtGrowthRateTotal;
        
        %%
        reactionTypeSum = sum(a);
        
        a0 = sum(reactionTypeSum);

        if a0 <= 0
            warning('No more events to execute');
            t = tend;
            break;
        end

        r1 = rand;
        r2 = rand;

        tau = (-log(r1)/a0);
        prop = a0 * r2;
        
        sum_j = 0;
        reactionFound = 0;

        for i = 1:nReactions
            inext = i;
            sum_j1 = sum_j;
            sum_j = sum_j + a(i);
            if sum_j >= prop && prop > sum_j1
                reactionFound = 1;
                break;
            end
        end
        
        if ~reactionFound
            error('Problem sampling reactions');
        end
        
        if t < warmupTime && t + tau > warmupTime
            t = warmupTime; % Would have gone over warm up time - stop and resample with new propensities
            continue;
        else
            t = t + tau;
        end
        
        delta = sum_j - prop;
        switch inext
            case 1
                % Next reaction is Kip2 attachment 
                site = find(attachmentCumSum >= delta, 1);
                if site > currentMtLength || isnan(mt(site))
                    error('Kip2 should not interact past MT end!');
                end
                if mt(site) ~= 0
                    error('Site already bound!');
                end
                mt(site) = 1;
                movingKip2s(site) = (site == currentMtLength) || (~mt(site+1));
                nKip2Free = nKip2Free - 1;
                if site > 1
                    movingKip2s(site - 1) = 0;
                end
            case 2
                % Next reaction is Kip2 stepping
                site = find(stepCumSum >= delta, 1);
                if site > currentMtLength || isnan(mt(site))
                    error('Kip2 should not interact past MT end!');
                end
                if mt(site) ~= 1
                    error('Kip2 location wrong!');
                end
                mt(site) = 0;
                movingKip2s(site) = 0;
                if site ~= currentMtLength
                    % Motor doesn't fall off
                    if mt(site + 1) ~= 0
                        error('Site already bound!');
                    end
                    mt( site + 1) = 1;
                    movingKip2s(site + 1) = ((site + 1) == currentMtLength) || (~mt(site + 2));
                else
                    % Motor falls off
                    %updateRuns(runLengthKip2(site) + 1, t + tau - attachedTimeKip2(site));
                    nKip2Free = nKip2Free + 1;
                end
                if site > 1
                    movingKip2s(site - 1) = mt(site - 1);
                end
            case 3
                % Next reaction is Kip2 unbinding in other zone
                site = find(detachmentCumSum >= delta, 1);
                if site > currentMtLength || isnan(mt(site))
                    error('Kip2 should not interact past MT end!');
                end
                if mt(site) ~= 1
                    error('Kip2 location wrong!');
                end

                mt(site) = 0;
                movingKip2s(site) = 0;
                %updateRuns(runLengthKip2(site), t + tau - attachedTimeKip2(site));
                if site > 1
                    movingKip2s(site - 1) = mt(site - 1);
                end
                nKip2Free = nKip2Free + 1;
            case 4
                % Next reaction elongates microtubule
                if mt(currentMtLength) == 1 && movingKip2s(currentMtLength) == 0
                    % Update Kip2 movement state (last Kip2 will now be able to move!)
                    movingKip2s(currentMtLength) = 1;
                end
                
                currentMtLength = currentMtLength + 1;
                mt(currentMtLength) = 0;
                movingKip2s(currentMtLength) = 0;
            case default
                error('Unknown reaction');
        end
        
        step = step + 1;
    end
    
    while tIndex <= length(tspan) && t > tspan(tIndex)
        mtState(tIndex, :) = mt;
        tState(tIndex) = t;
        kip2FreeState(tIndex) = nKip2Free;
        tIndex = tIndex + 1;
    end
    
    result = struct;
    result.mtState = mtState(1:tIndex-1, :);
    result.tState = tState(1:tIndex-1);
    result.steps = step;
    
    result.runtime = toc(localStart);
    result.parameters = parameters;
    result.tspan = tspan;
    
    parsave(currentFileName, result);
end
