function loadDependencies()
% loadDependencies: Downloads dependencies via Git and adds them to the MATLAB path
%
%   Usage:
%   loadDependencies()

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)


    if ~exist(['dependencies' filesep 'shadedErrorBar' filesep 'shadedErrorBar.m'], 'file')
        warning('Git submodules missing. Attempting to download them. Note that this requires installation of git on Linux/OS X or Git for Windows.');
        if ispc
            status = system('git submodule update --init --recursive');
        else
            status = system('LD_LIBRARY_PATH= && git submodule update --init --recursive');
        end
        if status ~= 0
            error('Could not initialize submodules');
        end
    end

    %%
    addpath('util');
    addpath(['dependencies' filesep 'interp1qr']);
    addpath(['dependencies' filesep 'violinPlot']);
    addpath(['dependencies' filesep 'shadedErrorBar']);
    addpath(['dependencies' filesep 'export_fig']);
    addpath(['dependencies' filesep 'legendflex' filesep 'legendflex']);
    addpath(['dependencies' filesep 'legendflex' filesep 'setgetpos_V1.2']);
    
end
