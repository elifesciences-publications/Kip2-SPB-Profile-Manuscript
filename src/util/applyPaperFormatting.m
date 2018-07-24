function applyPaperFormatting()
% applyPaperFormatting: Make everything Helvetica, font size 12

% © 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)
    fig = gcf;
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 12) 
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Helvetica')
end