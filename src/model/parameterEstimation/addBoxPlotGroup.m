function [x, g] = addBoxPlotGroup(x, g, x_new )
    if isempty(g)
        maxG = 1;
    else
        maxG = 1 + max(g);
    end
    x = [x; x_new];
    g = [g; zeros(size(x_new)) + maxG];
end

