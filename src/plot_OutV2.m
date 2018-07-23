function [] = plot_OutV2(OutV,paramSpecs,varargin)
    %% Set mode
    p = inputParser;
    addParameter(p, 'plotMode', 'sampleValues'); % other options: density, levelSets
    addParameter(p, 'bounds', 'fromParamSpecs'); % other options: fromData
    addParameter(p, 'costThreshold', Inf);
    addParameter(p, 'plotSamples', true);
    addParameter(p, 'nDoF', 0); 
    addParameter(p, 'nParameters', 0);
    addParameter(p, 'levelSets', [0.95, 0.75, 0.5, 0.25, 0.05]);
    addParameter(p, 'ticks', {});
    addParameter(p, 'showMaxLikelihoodPoint', true);
    addParameter(p, 'pointSize', 30);
    
    parse(p, varargin{:});



    %% plot cost

    if ~isfield(OutV,'colnames')
        OutV.colnames = paramSpecs.names;
    end
    NotProjectedIdx = ismember(paramSpecs.names,OutV.colnames);
    % paraoptestimated = paraopt(NotProjectedIdx); % to add fitted parameter point

    paramSpecsestimated = paramSpecs(NotProjectedIdx,:);
    data = OutV;
    if isfield(OutV,'V')
        data.rowmat = data.V;
    end
    paramNames = paramSpecsestimated.names;
    bmin = paramSpecsestimated.bmin;
    bmax = paramSpecsestimated.bmax;
    paramtoplot = paramNames;
    n = length(paramtoplot);

    indices = zeros(1,n);
    x = zeros(size(OutV.V, 1), n);
    xIndices = (OutV.cost <= p.Results.costThreshold);
    
    
    if p.Results.showMaxLikelihoodPoint
        [~, maxLikelihoodIndex] = min(OutV.cost);
    end
    
    xmin = zeros(1,n);
    xmax = zeros(1,n);
    for k = 1:n
        indices(k) = find(ismember(paramNames,paramtoplot{k}));
        x(:,k) = OutV.V(:,k);
        
        switch p.Results.bounds
            case 'fromParamSpecs'
                xmin(k) = bmin(indices(k));
                xmax(k) = bmax(indices(k));
            case 'fromData' 
                xmin(k) = min(x(xIndices,k));% bmin(indices(k));
                xmax(k) = max(x(xIndices,k)); %bmax(indices(k));
        end
    end

    f = figure();

    factor = 1.25;
    cMap = parula();
    switch p.Results.plotMode
        case 'levelSets'
            cMap = flipud(cMap); 
        case 'sampleValues'
            cMap = flipud(cMap); 
    end
    
    
    %% Diagonal Subplots
    minCost = min(OutV.cost);
    if strcmp(p.Results.plotMode, 'levelSets') || strcmp(p.Results.plotMode, 'sampleValues') 
        maxCostValues = nan(size(p.Results.levelSets));
        for j = 1:length(p.Results.levelSets)
            maxCostValues(j) = chi2inv(p.Results.levelSets(j), p.Results.nDoF - p.Results.nParameters);
        end
        maxCost = max(maxCostValues);
    else
        maxCost = max(OutV.cost(xIndices));
    end
    
    for diagonalIndex = 1:n
        subplot(n,n,(diagonalIndex-1)*(n+1) + 1);
        hold on;
        
        levelSetIndices = OutV.cost <= maxCost;
            
        scatter(x(levelSetIndices, diagonalIndex), OutV.cost(levelSetIndices), p.Results.pointSize, OutV.cost(levelSetIndices), 'filled' );
        

        if p.Results.showMaxLikelihoodPoint
            scatter(x(maxLikelihoodIndex, diagonalIndex), OutV.cost(maxLikelihoodIndex), p.Results.pointSize, OutV.cost(maxLikelihoodIndex), 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 1);
        end
        pos = get(gca, 'Position');
        pos(3) = 0.7*factor*pos(3);
        pos(4) = factor*pos(4);
        set(gca, 'Position', pos);      
        foo2 = gca;
        foo2.XLabel.Color = [0 0 0];
        foo2.YLabel.Color = [0 0 0];
        foo2.XColor = [0 0 0];
        foo2.YColor = [0 0 0];
        axis([xmin(diagonalIndex) xmax(diagonalIndex) minCost maxCost]);
        ylabel('Cost');
        box off

        if paramSpecs.islog(diagonalIndex)
            set(gca, 'XScale', 'log');
        end
        if diagonalIndex ~= n
            set(gca, 'XTick',[]);
        else
            if ~isempty(p.Results.ticks)
                set(gca, 'XTick', p.Results.ticks{diagonalIndex});
            end
            xlabel(paramSpecsestimated.legend{diagonalIndex});
        end
    end
    %% Off-diagonal subplots
    
    for row = 2:n
        for col = 1:row-1
            subplot(n,n,(row-1)*n + col);
            %
            xSpan = linspace(xmin(col), xmax(col), 101);
            ySpan = linspace(xmin(row), xmax(row), 101);
            hold on;
            switch p.Results.plotMode
                case 'sampleValues'
                    scatter(x(:, col), x(:, row), p.Results.pointSize, 'filled', 'MarkerFaceAlpha', 0.05, 'MarkerFaceColor', [0.5 0.5 0.5]);
                    
                    for j = 1:length(p.Results.levelSets)
                        levelSetThreshold = maxCostValues(j);
                        levelSetIndices = OutV.cost <= levelSetThreshold;
                        if j < length(p.Results.levelSets)
                            levelSetIndices = levelSetIndices & (OutV.cost > maxCostValues(j + 1));
                        end
                        
                        xSet = x(levelSetIndices, col);

                        
                        if ~isempty(xSet)
                            ySet = x(levelSetIndices, row);
                            costSet = OutV.cost(levelSetIndices);
                            normalizedCost = levelSetThreshold - minCost;
                            normalizedCost = normalizedCost ./ (maxCost - minCost);

                            scatter(xSet, ySet, p.Results.pointSize, costSet, 'filled', 'MarkerFaceAlpha', (1 - normalizedCost) * 0.5 + 0.5);
                        end
                    end

                case 'density'
                    [~,density,X,Y] = kde2d([x(xIndices,col),x(xIndices,row)],10);
                    if ~isempty(density)
                        contour(X,Y,density,'LineWidth',1);
                    end
                case 'levelSets'
                    if p.Results.plotSamples
%                         normalizedCost = OutV.cost(OutV.cost <= maxCost) - minCost;
%                         normalizedCost = normalizedCost ./ (maxCost - minCost);
%                         normalizedCost = max(normalizedCost, 0.01);
                        scatter(x(OutV.cost <= maxCost,col), x(OutV.cost <= maxCost,row), 10, OutV.cost(OutV.cost <= maxCost), 'filled');
                        %scatter(x(OutV.cost <= maxCost,col), x(OutV.cost <= maxCost,row), 30 * normalizedCost, OutV.cost(OutV.cost <= maxCost));
                    end
                    for j = 1:length(p.Results.levelSets)
                        levelSetThreshold = maxCostValues(j);
                        levelSetIndices = OutV.cost <= levelSetThreshold;
                        
                        xSet = x(levelSetIndices, col);
                        ySet = x(levelSetIndices, row);
                        costSet = OutV.cost(levelSetIndices);
                        if length(xSet) > 2
                            
                            
                            hullIndices = boundary(xSet, ySet, 1);

                            SVMModel = fitcsvm([x(:, col) x(:, row)], OutV.cost <= levelSetThreshold, 'Cost',[0 1;100 0], 'KernelFunction', 'gaussian')
                            [X,Y] = meshgrid(xSpan, ySpan);
                            label = predict(SVMModel,[X(:) Y(:)]);
                            any(label)
                            
                            label = double(reshape(label, 101, 101));
                            contour(X, Y, label, 'LineWidth', 1);
                            
                            colorIndex = 1 + round((levelSetThreshold - minCost) / (maxCost - minCost) * (size(cMap, 1) - 1));
                            plot(xSet(hullIndices), ySet(hullIndices), 'LineWidth', 1, 'Color', cMap(colorIndex, :));
                            
%                             xLowRes = [xSet(hullIndices) ySet(hullIndices)]';
% 
%                             indicesToUse = linspace(1, length(xLowRes),1001);
%                             indicesToUse = (indicesToUse >= 5) & (indicesToUse <= size(xLowRes,2) - 4);
% 
%                             xHighRes = pchip(linspace(1,length(xLowRes),size(xLowRes,2)), xLowRes, linspace(1, size(xLowRes,2),1001));
%                             %contour(X, Y, density, 'LineWidth', 1);
% 
%                             colorIndex = round((levelSetThreshold - minCost) / (maxCost - minCost) * size(cMap, 1));
%                             hold on;
%                             plot(xHighRes(1,indicesToUse), xHighRes(2,indicesToUse), 'LineWidth', 1, 'Color', cMap(colorIndex, :));
                        end
                    end

                    
                    %plot(
                    %contour
            end
            
            if p.Results.showMaxLikelihoodPoint
                scatter(x(maxLikelihoodIndex, col), x(maxLikelihoodIndex, row), p.Results.pointSize, OutV.cost(maxLikelihoodIndex), 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 1);
                plot([xmin(col) x(maxLikelihoodIndex, col)], x(maxLikelihoodIndex, row) * [1 1], '--r');
                plot(x(maxLikelihoodIndex, col) * [1 1], [xmin(row) x(maxLikelihoodIndex, row)], '--r');
                %plot(x(maxLikelihoodIndex, col), x(maxLikelihoodIndex, row), '.', 'Color', 'red', 'MarkerSize', 15);
            end
           

            pos = get(gca, 'Position');
            pos(3) = 0.7*factor*pos(3);
            pos(4) = factor*pos(4);
            set(gca, 'Position', pos);
            foo2 = gca;
            foo2.XLabel.Color = [0 0 0];
            foo2.YLabel.Color = [0 0 0];
            foo2.XColor = [0 0 0];
            foo2.YColor = [0 0 0];            
            axis([xmin(col) xmax(col) xmin(row) xmax(row)]);
            
            if col ~= 1
                set(gca, 'YTick',[]); 
            else
                ylabel(paramSpecsestimated.legend{row});
                if ~isempty(p.Results.ticks)
                    set(gca, 'YTick', p.Results.ticks{row});
                end
            end
            if row ~= n
                set(gca, 'XTick',[]);
            else
                if ~isempty(p.Results.ticks)
                    set(gca, 'XTick', p.Results.ticks{col});
                end
                xlabel(paramSpecsestimated.legend{col});
            end
            
            if paramSpecs.islog(col)
                set(gca, 'XScale', 'log');
            end
            
            if paramSpecs.islog(row)
                set(gca, 'YScale', 'log');
            end
        end

    end
    colormap(cMap);
    
    c = colorbar();
    c.XTick = fliplr(maxCostValues);
    c.Limits(2) = maxCost;
    c.Color = [0 0 0];
    for i = 1:length(maxCostValues)
        c.TickLabels{i} = sprintf('%g%%', 100*p.Results.levelSets(length(maxCostValues) - i + 1));
    end

    pos = get(c, 'Position');
    
    
    set (c, 'Position', [pos(1)+0.75/n pos(2)+1/n pos(3) pos(4)*(n-1)]);
    
    switch p.Results.plotMode
        case 'density'
            c.Label.String = 'Probability';
        case 'sampleValues'
            c.Label.String = 'Confidence Region';
        case 'levelSets'
            c.Label.String = 'Confidence Region';
    end

    if nargin == 3
        set(f,'PaperOrientation','landscape','PaperUnits','Normalized','PaperPosition',[0 0 1 1]);
        print(fileName, '-dpdf');
    end
end