%function determinePlusEndNumbers
    myFolder = fileparts(mfilename('fullpath'));

    figureResultFolder = [myFolder filesep '..' filesep '..' filesep 'figures' filesep 'model']; 
    rng('default');
    %%
    tspan = (7200+(0:10:14400));


    maxLength       = 268; % binding sites on protofilament
    nKip2Free       = 80; % Kip2 molecules per protofilament
    k_on            = 6.1054e-04; % s^-1 nM Kip2 ^ -1
    k_in            = 0.3162; % s^-1 nM Kip2 ^ -1
    k_out           = 3.0066; % % 1/s

    v_step_Kip2only = 6.259 * 1000 / 60; %3.5 * 1000 / 60; % um/min -> nm/s, Xiuzhen's data
    k_step_Kip2only = v_step_Kip2only / 8; % nm/s -> sites/s

    currentParameters                  = struct;
    currentParameters.maxLength        = maxLength;
    currentParameters.nKip2Free        = nKip2Free;
    currentParameters.k_on_mt          = ones(1, currentParameters.maxLength) * k_on;
    currentParameters.k_on_mt(1)       = k_on + k_in;
    currentParameters.k_step_mt        = ones(1, currentParameters.maxLength) * k_step_Kip2only;
    currentParameters.k_step_mt(end)   = 0;
    currentParameters.k_detach_mt      = zeros(1, currentParameters.maxLength);
    currentParameters.k_detach_mt(end) = k_out;
    currentParameters.reportRunLengths = false;
    %%
    parameterIndex = 0;
    currentFileName = '';
    currentOutFileName = '';

    simResult   = simulateModel(currentParameters, parameterIndex, tspan, currentFileName);

    vmResult = computeMeasurementModel(simResult, currentOutFileName);
    result = vmResult;
    %%
    meanBindingEnd = mean(simResult.mtState(:,end));
    fprintf('Mean Kip2 molecules at last site (for 13 protofilaments): %g\n', round(meanBindingEnd * 13));
    fprintf('Mean Kip2 molecules on the last 500 nm (for 13 protofilaments): %g\n', round(mean(sum(simResult.mtState(:,end-62:end),2))*13, -1)); % last 500 nm -> 63 sites
    
    %foo = mean(simResult.mtState,1);
    %figure;
    %plot(foo);
    return;


