function plotViolinPlots(prefix, parameterCombinations, pMin, pMax, exportFigures, parameterEstimateFigureResultFolder)
    %% Plot beeswarm plots
    figure;
    currAx = gca();
    colors = currAx.ColorOrder;
    [xConc, gConc] = addBoxPlotGroup([], [], parameterCombinations(1, :) ./ 140 .* 61);
    violins = violinplot(xConc, gConc, 'Bandwidth', 30, 'Width', 0.4, 'ShowData', false);
    hold on;
    xlim([0.5 1.5]);
    plot([0.6 1.4], [1 1] .* pMin(1) ./ 140 .* 61, '--k');
    plot([0.6 1.4], [1 1] .* pMax(1) ./ 140 .* 61, '--k');
    ylabel('nM');
    ylim([0 230]);
    
    currAx.XTickLabel = {'[Kip2]_{total}'};
    currAx.XLabel.Color = [0 0 0];
    currAx.YLabel.Color = [0 0 0];
    currAx.XColor = [0 0 0];
    currAx.YColor = [0 0 0];

    applyPaperFormatting

    currFig = gcf();
    currFig.Position(3) = 0.25*currFig.Position(3);
    currFig.Color = 'white';
    
    concentrations = (parameterCombinations(1, :) ./ 140 .* 61);
    fprintf('Mean [Kip2]: %g,\tStd:%g\n', mean(concentrations), std(concentrations));
    fprintf('Median [Kip2]: %g\n', median(concentrations));

    if exportFigures
        print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_conc.pdf'], '-dpdf');
        export_fig([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_conc.png'], '-q101', '-r300');
    end


    figure;
    xConc = [];
    gConc = [];
    [xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(parameterCombinations(2, :)));
    [xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(parameterCombinations(3, :)));
    violins = violinplot(xConc, gConc, 'Bandwidth', 1, 'Width', 0.4, 'ShowData', false);
    hold on;
    plot([0.6 1.4], log10([1 1] .* pMin(2)), '--k');
    plot([0.6 1.4], log10([1 1] .* pMax(2)), '--k');

    plot([1.6 2.4], log10([1 1] .* pMin(3)), '--k');
    plot([1.6 2.4], log10([1 1] .* pMax(3)), '--k');
    
    %violins(2).ViolinColor = colors(2, :);
    xlim([0.5 2.5]);
    ylim([-5.5 2]);
    ylabel('nM Kip2^{-1} s^{-1}');

    currAx = gca();
    currAx.XTickLabel = {'k_{on}', 'k_{in}'};
    currAx.YTick = -5:2;
    currAx.YTickLabel = {'10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{-0}', '10^{1}', '10^{2}'};
    currAx.XLabel.Color = [0 0 0];
    currAx.YLabel.Color = [0 0 0];
    currAx.XColor = [0 0 0];
    currAx.YColor = [0 0 0];

    applyPaperFormatting

    currFig = gcf();
    currFig.Position(3) = 0.4*currFig.Position(3);
    currFig.Color = 'white';
    
    fprintf('Mean k_on: %g,\t%g\n', mean(parameterCombinations(2, :)), std(parameterCombinations(2, :)));
    fprintf('Median k_on: %g\n', median(parameterCombinations(2, :)));
    
    fprintf('Mean k_in: %g\t%g\n', mean(parameterCombinations(3, :)), std(parameterCombinations(3, :)));
    fprintf('Median k_in: %g\n', median(parameterCombinations(3, :)));
    
    if exportFigures
        print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_kon_kin.pdf'], '-dpdf');
        export_fig([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_kon_kin.png'], '-q101', '-r300');
    end

    figure;
    xConc = [];
    gConc = [];
    [xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(parameterCombinations(2, :) .* parameterCombinations(1, :) ./ 140 .* 61));
    [xConc, gConc] = addBoxPlotGroup(xConc, gConc, log10(parameterCombinations(3, :) .* parameterCombinations(1, :) ./ 140 .* 61));
    violins = violinplot(xConc, gConc, 'Bandwidth', 1, 'Width', 0.4, 'ShowData', false);
    
    hold on;
    plot([0.6 1.4], log10([1 1] .* pMin(2) .* pMin(1) ./ 140 .* 61), '--k');
    plot([0.6 1.4], log10([1 1] .* pMax(2) .* pMax(1) ./ 140 .* 61), '--k');

    plot([1.6 2.4], log10([1 1] .* pMin(3) .* pMin(1) ./ 140 .* 61), '--k');
    plot([1.6 2.4], log10([1 1] .* pMax(3) .* pMax(1) ./ 140 .* 61), '--k');
    
    %violins(2).ViolinColor = colors(2, :);
    xlim([0.5 2.5]);
    ylim([-4 5]);
    ylabel('s^{-1}');

    currAx = gca();
    currAx.XTickLabel = {'r_{on, max}', 'r_{in, max}'};
    currAx.YTick = -5:5;
    currAx.YTickLabel = {'10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{-0}', '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}'};
    currAx.XLabel.Color = [0 0 0];
    currAx.YLabel.Color = [0 0 0];
    currAx.XColor = [0 0 0];
    currAx.YColor = [0 0 0];

    applyPaperFormatting

    currFig = gcf();
    currFig.Position(3) = 0.4*currFig.Position(3);
    currFig.Color = 'white';
    
    fprintf('Mean r_on,max: %g,\tStd: %g\n', mean(parameterCombinations(1, :) .* parameterCombinations(2, :)), std(parameterCombinations(1, :) .* parameterCombinations(2, :)));
    fprintf('Median r_on,max: %g\n', median(parameterCombinations(1, :) .* parameterCombinations(2, :)));
    
    fprintf('Mean r_in,max: %g,\tStd: %g\n', mean(parameterCombinations(1, :) .* parameterCombinations(3, :)), std(parameterCombinations(1, :) .* parameterCombinations(3, :)));
    fprintf('Median r_in,max: %g\n', median(parameterCombinations(1, :) .* parameterCombinations(3, :)));

    if exportFigures
        print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_ron_rin.pdf'], '-dpdf');
        export_fig([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_ron_rin.png'], '-q101', '-r300');
    end

    figure;
    [xConc, gConc] = addBoxPlotGroup([], [], parameterCombinations(4, :));
    violinplot(xConc, gConc, 'Bandwidth', 0.5, 'Width', 0.4, 'ShowData', false);
    xlim([0.5 1.5]);
    ylim([0.5 14]);
    ylabel('s^{-1}');


    currAx = gca();
    currAx.XTickLabel = {'k_{out}'};
    currAx.XLabel.Color = [0 0 0];
    currAx.YLabel.Color = [0 0 0];
    currAx.XColor = [0 0 0];
    currAx.YColor = [0 0 0];
    
    plot([0.6 1.4], [1 1] .* pMin(4), '--k');
    plot([0.6 1.4], [1 1] .* pMax(4), '--k');
    
    applyPaperFormatting

    currFig = gcf();
    currFig.Position(3) = 0.25*currFig.Position(3);
    currFig.Color = 'white';
    
    fprintf('Mean k_out: %g\n', mean(parameterCombinations(4, :)));
    fprintf('Median k_out: %g\n', median(parameterCombinations(4, :)));    

    if exportFigures
        print([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_kout.pdf'], '-dpdf');
        export_fig([parameterEstimateFigureResultFolder filesep 'pdf' filesep prefix '_violin_kout.png'], '-q101', '-r300');
    end
end