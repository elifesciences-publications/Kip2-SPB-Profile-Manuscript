function exportTwoColorCondtionResultToData(conditionResult, fileName)
    dataTableGreen = array2table([conditionResult.lengthVector conditionResult.greenIntensities]);
    dataTableGreen.Properties.VariableNames{1} = 'x';
    dataTableGreen.Properties.VariableNames(2:end) = conditionResult.cellNames;
    writetable(dataTableGreen, [fileName '-green.txt'], 'Delimiter', '\t')
    
    dataTableRed = array2table([conditionResult.lengthVector conditionResult.redIntensities]);
    dataTableRed.Properties.VariableNames{1} = 'x';
    dataTableRed.Properties.VariableNames(2:end) = conditionResult.cellNames;
    writetable(dataTableRed, [fileName '-red.txt'], 'Delimiter', '\t')
end

