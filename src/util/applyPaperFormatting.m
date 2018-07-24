function applyPaperFormatting()
    fig = gcf;
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 12) 
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica')
end