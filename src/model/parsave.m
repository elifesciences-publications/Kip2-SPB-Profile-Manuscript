function parsave(fname, result)
    if isempty(fname)
        %warning('No filename specified - not saving file');
        return
    end
    save(fname, 'result', '-v7.3');
end