function [fitresult, gof] = laserModelFit(quantilesAllLow, quantilesAllHigh, doPlot)


%  Create a fit.
%
%  Data for fit:
%      X Input : quantilesAllLow
%      Y Output: quantilesAllHigh
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Nov-2018 14:07:20


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( quantilesAllLow, quantilesAllHigh );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
%opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if doPlot 
    % Create a figure for the plots.
    figure( 'Name', 'untitled fit 1' );

    % Plot fit with data.
    subplot( 2, 1, 1 );
    h = plot( fitresult, xData, yData );
    legend( h, 'quantilesAllHigh vs. quantilesAllLow', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel quantilesAllLow
    ylabel quantilesAllHigh
    grid on

    % Plot residuals.
    subplot( 2, 1, 2 );
    h = plot( fitresult, xData, yData, 'residuals' );
    legend( h, 'untitled fit 1 - residuals', 'Zero Line', 'Location', 'NorthEast' );
    % Label axes
    xlabel quantilesAllLow
    ylabel quantilesAllHigh
    grid on
    
end

