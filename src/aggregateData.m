function aggregateData()
% aggregateData: Aggregates 5 px wide line profiles from individual green / red 
% channel CSV files, imports them into tables, and writes the tables into CSV
% files again.
%
%   Usage:
%   aggregateData()


% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)
baseDirectory = ['..' filesep 'rawData' filesep];
resultDirectory = ['..' filesep 'analyzedData' filesep 'aggregatedProfiles' filesep];

if ~exist(resultDirectory, 'dir')
    mkdir(resultDirectory);
end

directories = {
    '20180130-15100b-Kip2-3sfGFP-Spc42-mCherry-sum-projection'
    '20180130-15100-Kip2-3sfGFP-Spc42-mCherry-sum-projection'
    '20180212_14590b_Kip3-3sfGFP_Spc42-mCherry_sum_projection'    
    '20180212_14590_Kip3-3sfGFP_Spc42-mCherry_sum_projection'  
    '20180224_15100_Kip2-3sfGFP_Spc42-mCherry_sumProjecion'
    '20180224_15102_Kip2-S63A-3sfGFP_Spc42-mCherry_sumProjecion'
    };

outFileNames = {
    '20180130_15100b_Kip2-3sfGFP_Spc42-mCherry'
    '20180130_15100_Kip2-3sfGFP_Spc42-mCherry'
    '20180212_14590b_Kip3-3sfGFP_Spc42-mCherry'    
    '20180212_14590_Kip3-3sfGFP_Spc42-mCherry'
    '20180224_15100_Kip2-3sfGFP_Spc42-mCherry'
    '20180224_15102_Kip2-S63A-3sfGFP_Spc42-mCherry'
    };

for dirIndex = 1:length(directories)
    directory = [baseDirectory directories{dirIndex}];
    allGreenCellFiles = dir([directory filesep '**' filesep '*g.csv']);
    allRedCellFiles = dir([directory filesep '**' filesep '*r.csv']);
    
    if length(allGreenCellFiles) ~= length(allRedCellFiles)
        f1 = {allGreenCellFiles.name}';
        f2 = {allRedCellFiles.name}';
        
        maxInd = max(size(f1, 1), size(f2, 1));
        for i=1:maxInd
            if ~strcmp(strrep(f1{i},'g','r'), f2{i})
                 error('# green profiles should be = # red profiles, double-check quantificiation!\nFile "%s"', [allGreenCellFiles(i).folder filesep allGreenCellFiles(i).name]); 
            end
        end
		
		
        error('# green profiles should be = # red profiles, double-check quantificiation!');
    end

    greenFilesByCondition = struct;
    for i = 1:length(allGreenCellFiles)
        currentGreenFile = allGreenCellFiles(i);
        currentGreenFileName = [currentGreenFile.folder filesep currentGreenFile.name];
        
        if contains(currentGreenFileName, 'G1')
            if isfield(greenFilesByCondition, 'G1')
                greenFilesByCondition.G1(end+1) = allGreenCellFiles(i);
            else
                greenFilesByCondition.G1 = allGreenCellFiles(i);
            end    
        elseif contains(currentGreenFileName, 'Distal')
            if isfield(greenFilesByCondition, 'Distal')
                greenFilesByCondition.Distal(end+1) = allGreenCellFiles(i);
            else
                greenFilesByCondition.Distal = allGreenCellFiles(i);
            end           
        else
            if isfield(greenFilesByCondition, 'metaPhase')
                greenFilesByCondition.metaPhase(end+1) = allGreenCellFiles(i);
            else
                greenFilesByCondition.metaPhase = allGreenCellFiles(i);
            end
        end
        
    end

    %%
    
    cellCycles = fieldnames(greenFilesByCondition)';
    
    for cellCycle = cellCycles
        T_green = table([0],'VariableNames',{'X'});
        T_red   = table([0],'VariableNames',{'X'});

        currentGreenFilesInCellCycle = greenFilesByCondition.(cellCycle{1});
        for i = 1:length(currentGreenFilesInCellCycle)
            currentGreenFile = currentGreenFilesInCellCycle(i);

            currentGreenFileName = [currentGreenFile.folder filesep currentGreenFile.name];
            currentRedFileName = currentGreenFileName;
            currentRedFileName(end-4) = 'r';

            if ~exist(currentRedFileName, 'file')
                error([currentRedFileName ' does not exist!']);
            end
            if contains(currentGreenFileName, 'DONOTUSE', 'IgnoreCase', true)
                warning(['Not using file "' currentFileName '"']);
                continue;
            end
            
            [temp, dirName] = fileparts(currentGreenFile.folder);
            
            if strcmpi(cellCycle, dirName)
            	[~, dirName] = fileparts(temp); % strip G1 / Distal
            end
            varName = ['Y_' dirName '_' currentGreenFile.name(1:end-4)];
            varName = strrep(varName, '-', '_');
            
            T_g = readtable(currentGreenFileName);
            T_g.Properties.VariableNames{1} = ['X'];
            T_g.X = round(T_g.X, 4);
            
            T_g.Properties.VariableNames{2} = varName;
            
            for j = 1:min(size(T_green.X,1), size(T_g.X,1))
                if abs(T_green.X(j) - T_g.X(j) < 0.05)
                    T_g.X(j) = T_green.X(j);
                end
             end
            T_green = outerjoin(T_green, T_g, 'Keys','X', 'MergeKeys',true);
            
            T_r = readtable(currentRedFileName);
            T_r.Properties.VariableNames{1} = ['X'];
            T_r.X = round(T_r.X, 4);
            for j = 1:min(size(T_red.X,1), size(T_r.X,1))
                if abs(T_red.X(j) - T_r.X(j) < 0.05)
                    T_r.X(j) = T_red.X(j);
                end
             end
            T_r.Properties.VariableNames{2} = varName;
            T_red = outerjoin(T_red, T_r, 'Keys','X', 'MergeKeys',true);

            if any(size(T_g) ~= size(T_r))
                error(['Size mismatch: "' currentGreenFileName '"']);
            end

            if all(all(T_r{:,:} == T_g{:,:}))
                error(['Red and green profile are the same: ' currentGreenFileName]);
            end
        end

        writetable(T_green, [resultDirectory outFileNames{dirIndex} '-' cellCycle{1} '-named-green.txt'], 'Delimiter', '\t')
        writetable(T_red,   [resultDirectory outFileNames{dirIndex} '-' cellCycle{1} '-named-red.txt'], 'Delimiter', '\t')
    end
end