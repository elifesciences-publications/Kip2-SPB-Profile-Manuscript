function parsave(fileName, result)
    if isempty(fileName)
        return
    end
    save(fileName, 'result', '-v7.3');
end