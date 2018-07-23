function [conditions, maxLength] = displayConditions(conditionResults)
    conditions = {};

    maxLength = 0;
    for currentConditionIndex = 1:length(conditionResults)
        c = conditionResults{currentConditionIndex};
        conditions{end+1} = c.condition;
        meanLength = mean(c.lengthsUM); 
        if meanLength > maxLength
            maxLength = meanLength;
        end
        fprintf('%2i: %s (n = %i)\n', currentConditionIndex, c.condition, size(c.lengths, 1));
    end
end