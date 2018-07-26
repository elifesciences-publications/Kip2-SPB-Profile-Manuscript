function parsave(fileName, result)
% parsave: wrap save function such that it can be used inside a parfor loop.
% 
%   Usage:
%   parsave(fileName, result)
%
% Arguments:
%   fileName: char array
%     fileName to save to
%   result: var
%     variable to save

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

    if isempty(fileName)
        return
    end
    save(fileName, 'result', '-v7.3');
end