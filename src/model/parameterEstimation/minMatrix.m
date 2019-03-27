function [minimum, rows, cols] = minMatrix(mat)
    [minimum, minimumIndex] = min(mat(:));
    
    [rows, cols] = ind2sub(size(mat), minimumIndex);
end