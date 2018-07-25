function mergeAndScaleDatasetsTwoColor(file1, scaling1, file2, scaling2, outFile)
% mergeAndScaleDatasetsTwoColor: Merge red/green datasets while applying 
% an affine scaling, and export to outFile.green / outFile.red CSV files
% 
%   Usage:
%   mergeAndScaleDatasetsTwoColor(file1, scaling1, file2, scaling2, outFile)
%
% Arguments:
%   file1/file2: struct
%     .green: input green fluorescence CSV file
%     .red  : input red fluorescence CSV file
%   scaling1/2: struct 
%     .offset: offset applied in the affine scaling
%     .factor: scaling factor for the affine scaling
%   outFile: struct
%     .green: output green fluorescence CSV file
%     .red  : output red fluorescence CSV file

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)
    
    for channel = {'green', 'red'}
        currentChannel = channel{1};
        T1.(currentChannel) = readtable(file1.(currentChannel));
        T2.(currentChannel) = readtable(file2.(currentChannel));

        T1.(currentChannel).X = round(T1.(currentChannel).X, 4);
        for j = 1:min(size(T1.(currentChannel).X, 1), size(T2.(currentChannel).X, 1))
            if abs(T1.(currentChannel).X(j) - T2.(currentChannel).X(j) < 0.05)
                T2.(currentChannel).X(j) = T2.(currentChannel).X(j);
            end
        end

        T1.(currentChannel){:, 2:end} = scaling1.(currentChannel).offset + scaling1.(currentChannel).factor .* T1.(currentChannel){:, 2:end};
        T2.(currentChannel){:, 2:end} = scaling2.(currentChannel).offset + scaling2.(currentChannel).factor .* T2.(currentChannel){:, 2:end};

        T.(currentChannel) = outerjoin(T1.(currentChannel), T2.(currentChannel), 'Keys', 'X', 'MergeKeys', true);
        writetable(T.(currentChannel), outFile.(currentChannel), 'Delimiter', '\t');
    end
end