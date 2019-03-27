function plotEstimateErrorBar(parameterCombinations, uniqueCompatibleParameterSets, p, parameterIndex, x, color)
    plot(x, p(parameterIndex), 'x', 'Color', color, 'LineWidth', 0.8);
    if ~isempty(uniqueCompatibleParameterSets)
        uniqueParameterValues = unique(parameterCombinations(parameterIndex, uniqueCompatibleParameterSets));
        y = p(parameterIndex); 
        neg = y - uniqueParameterValues(1);
        pos = uniqueParameterValues(end) - y;
        errorbar(x,y,neg,pos,'Color', color, 'LineWidth', 0.8);
    end
end