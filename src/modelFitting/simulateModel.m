function result = simulateModel(parameters, parameterIndex, tspan, currentFileName)
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
    maxLength         = parameters.maxLength;
    
    k_on_mt           = parameters.k_on_mt;
    k_step_mt         = parameters.k_step_mt;
    k_detach_mt       = parameters.k_detach_mt;
    
    nKip2Free         = parameters.nKip2Free;
    
    reportRunLengths = parameters.reportRunLengths;
    
    
    tend = tspan(end); % seconds
    mt = zeros(1, maxLength);
    
    
    runLengthKip2 = zeros(1, maxLength);
    attachedTimeKip2 = zeros(1, maxLength);

    nReactions = 3;
    a = zeros(1, nReactions);
    
    step = 0;
    
    t = 0;
    tIndex = 1;
    mtState = zeros(length(tspan), maxLength);

    
    runSize = 100;
    nRuns = 0;

    runLengthsKip2 = nan(1, runSize);
    runTimesKip2   = nan(1, runSize);
    
    function updateRuns(runLength, runTime)
        if ~reportRunLengths
            return
        end
        nRuns = nRuns + 1;
        
        if nRuns > runSize
            runLengthsKip2 = [runLengthsKip2 nan(1, runSize)];
            runTimesKip2   = [runTimesKip2 nan(1, runSize)];
            runSize = runSize + runSize;
        end
        
        runLengthsKip2(nRuns) = runLength;
        runTimesKip2(nRuns) = runTime;
    end
    
    tState = zeros(size(tspan));
    while t < tend
        
        if nKip2Free < 0 || nKip2Free > parameters.nKip2Free
            error('Molecules not conserved!');
        end
        
        % 0: Free
        % 1: Kip2
        
        while tIndex <= length(tspan) && t > tspan(tIndex)
            mtState(tIndex, :) = mt;
            tState(tIndex) = t;
            tIndex = tIndex + 1;
        end

        %% Attachment       
        freeSites = mt == 0;
        
        r_on_mt = k_on_mt * (60/140) * nKip2Free; % Conversion: 140 molecules =^= 60 nM
        attachmentRates = double(freeSites) .* r_on_mt;
        attachmentCumSum = cumsum(attachmentRates);
        
        a(:, 1) = attachmentCumSum(end);
        
        %% Stepping
        kip2Motors = mt == 1;
        
        movingKip2s = kip2Motors & [freeSites(2:end) 1];
        
        % Motors 
        stepRates = double(movingKip2s) .* k_step_mt;
        stepCumSum = cumsum(stepRates);
        
        a(:, 2) = stepCumSum(end);
        
        %% Detachment   
        detachmentRates = double(kip2Motors) .* k_detach_mt;
        detachmentCumSum = cumsum(detachmentRates);
        a(:, 3) = detachmentCumSum(end);     
        
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
        jnext = 1;
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
        
        delta = sum_j - prop;
        switch inext
            case 1
                % Next reaction is Kip2 attachment 
                site = find(attachmentCumSum >= delta, 1);
                if mt(jnext, site) ~= 0
                    error('Site already bound!');
                end
                mt(jnext, site) = 1;
                runLengthKip2(jnext, site) = 0;
                attachedTimeKip2(jnext, site) = t + tau;
                nKip2Free = nKip2Free - 1;
            case 2
                % Next reaction is Kip2 stepping
                site = find(stepCumSum >= delta, 1);
                if mt(jnext, site) ~= 1
                    error('Kip2 location wrong!');
                end
                
                mt(jnext, site) = 0;
                if site ~= maxLength
                    % Motor doesn't fall off
                    if mt(jnext, site + 1) ~= 0
                        error('Site already bound!');
                    end
                    mt(jnext, site + 1) = 1;
                    runLengthKip2(jnext, site + 1) = runLengthKip2(jnext, site) + 1; 
                    attachedTimeKip2(jnext, site + 1) = attachedTimeKip2(jnext, site);
                else
                    % Motor falls off
                    updateRuns(runLengthKip2(jnext, site) + 1, t + tau - attachedTimeKip2(jnext, site));
                    nKip2Free = nKip2Free + 1;
                end
                runLengthKip2(jnext, site) = 0;
                attachedTimeKip2(jnext, site) = 0;
            case 3
                % Next reaction is Kip2 unbinding in other zone
                site = find(detachmentCumSum >= delta, 1);

                if mt(jnext, site) ~= 1
                    error('Kip2 location wrong!');
                end
                mt(jnext, site) = 0;
                
                updateRuns(runLengthKip2(jnext, site), t + tau - attachedTimeKip2(jnext, site));
                
                nKip2Free = nKip2Free + 1;
                attachedTimeKip2(jnext, site) = 0;
                runLengthKip2(jnext, site) = 0;
            case default
                error('Unknown reaction');
        end
        
        
        t = t + tau;
        
        step = step + 1;
    end
    
    while tIndex <= length(tspan) && t > tspan(tIndex)
        mtState(tIndex, :) = mt;
        tState(tIndex) = t;
        tIndex = tIndex + 1;
    end
    
    result = struct;
    result.mtState = mtState;
    result.runLengthsKip2 = runLengthsKip2(1:nRuns);
    result.runTimesKip2 = runTimesKip2(1:nRuns);
    result.steps = step;
    
    result.runtime = toc(localStart);
    result.parameters = parameters;
    result.tspan = tspan;
    result.parameterIndex = parameterIndex;
    
    parsave(currentFileName, result);
end
