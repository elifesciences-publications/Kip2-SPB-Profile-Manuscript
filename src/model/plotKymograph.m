function plotKymograph
    myFolder = fileparts(mfilename('fullpath'));

    figureResultFolder = [myFolder filesep '..' filesep '..' filesep 'figures' filesep 'model']; 
    rng('default');
    %%
    tspan = (7200+(0:1:87));


    maxLength       = 283; % binding sites on protofilament
    nKip2Free       = 80; % Kip2 molecules per protofilament
    k_on            = 6.1054e-04; % s^-1 nM Kip2 ^ -1
    k_in            = 0.3162; % s^-1 nM Kip2 ^ -1
    k_out           = 3.0066; % % 1/s

    v_step_Kip2only = 6.259 * 1000 / 60; %3.5 * 1000 / 60; % um/min -> nm/s, Xiuzhen's data
    k_step_Kip2only = v_step_Kip2only / 8; % nm/s -> sites/s

    currentParameters                  = struct;
    currentParameters.maxLength        = maxLength;
    currentParameters.nKip2Free        = nKip2Free;
    currentParameters.k_on_mt          = ones(1, currentParameters.maxLength) * k_on;
    currentParameters.k_on_mt(1)       = k_on + k_in;
    currentParameters.k_step_mt        = ones(1, currentParameters.maxLength) * k_step_Kip2only;
    currentParameters.k_step_mt(end)   = 0;
    currentParameters.k_detach_mt      = zeros(1, currentParameters.maxLength);
    currentParameters.k_detach_mt(end) = k_out;
    currentParameters.reportRunLengths = false;
    %%
    parameterIndex = 0;
    currentFileName = '';
    currentOutFileName = '';

    simResult   = simulateModel(currentParameters, parameterIndex, tspan, currentFileName);

    vmResult = computeMeasurementModel(simResult, currentOutFileName);
    result = vmResult;
    %%

    figure();

    j = 1;
    imnew = result.simResult.mtState'; %imresize(result.vm.binnedProfiles{j}(:,:)', 4);
    image(tspan, ((1:size(result.simResult.mtState,2))-1) .* 0.008,imnew, 'CDataMapping', 'scaled');
    hold on;
    currentAx = gca;
    currentAx.YDir = 'normal';
    colormap('gray');
    currentAx.YTick = [0 maxLength*0.008];
    currentAx.YTickLabel = {'SPB', '+ end'};
    currentAx.XTick = [];
    currentAx.TickDir = 'out';
    currentAx.Color = 'black';
    currentFig = gcf;
    currentFig.Position(4) = 0.5 * currentFig.Position(4) / 1.4;
    currentFig.Position(3) = currentFig.Position(3) * 1.1;

    plot([tspan(4) tspan(4)],[result.vm.binnedXCoords{j}(end-1) (result.vm.binnedXCoords{j}(end-1)-2) ],'w-','linewidth',2); 
    plot([tspan(7) tspan(7)+10],[result.vm.binnedXCoords{j}(end-1) result.vm.binnedXCoords{j}(end-1)],'w-','linewidth',2); 

    spaceText = text(tspan(9),result.vm.binnedXCoords{j}(end-1) - 1,'2 µm','horiz','center','vert','top'); 
    spaceText.Color = 'white';
    spaceText.Rotation = -90;


    timeText = text(tspan(7)+5,result.vm.binnedXCoords{j}(end-1)-0.05,'10 s','horiz','center','vert','top'); 
    timeText.Color = 'white';

    ylim([min(result.vm.binnedXCoords{j}) max(result.vm.binnedXCoords{j})]);
    applyPaperFormatting

    %%

    figure();

    j = 1;
    imnew = result.vm.binnedProfiles{j}(:,:)'; %imresize(result.vm.binnedProfiles{j}(:,:)', 4);
    image(tspan, result.vm.binnedXCoords{j} ,imnew, 'CDataMapping', 'scaled');
    hold on;
    currentAx = gca;
    currentAx.YDir = 'normal';
    colormap('gray');
    currentAx.YTick = [0 maxLength*0.008];
    currentAx.YTickLabel = {'SPB', '+ end'};
    currentAx.XTick = [];
    currentAx.TickDir = 'out';
    currentFig = gcf;
    currentFig.Position(4) = 0.5 * currentFig.Position(4) / 1.2;
    currentFig.Position(3) = currentFig.Position(3) * 1.1;

    axWidthRel = currentAx.Position(3) - currentAx.Position(1);
    axHeightRel = currentAx.Position(4) - currentAx.Position(2);

    yUm = max(result.vm.binnedXCoords{j}) - min(result.vm.binnedXCoords{j});

    figureRatio = 80.7/29.2;
    figureRatio = figureRatio * (3.9 / yUm);

    currentFig.Position(3) = figureRatio/ (axWidthRel / (currentFig.Position(4) * axHeightRel));
    currentFig.Position(3) * axWidthRel / (currentFig.Position(4) * axHeightRel);

    plot([tspan(2) tspan(2)],[result.vm.binnedXCoords{j}(end-1) (result.vm.binnedXCoords{j}(end-1)-2) ]-0.5,'Color', [244 234 14]./255,'linewidth',2.2); 
    plot([tspan(5) tspan(5)+10],[result.vm.binnedXCoords{j}(end-1) result.vm.binnedXCoords{j}(end-1)]-2.7,'Color', [244 234 14]./255,'linewidth',2.2); 

    spaceText = text(tspan(6),result.vm.binnedXCoords{j}(end-1) - 1.5,'2 µm','horiz','center','vert','top'); 
    spaceText.Color = [244 234 14]./255;
    spaceText.Rotation = -90;


    timeText = text(tspan(5)+5,result.vm.binnedXCoords{j}(end-1)-2.2,'10 s','horiz','center','vert','top'); 
    timeText.Color = [244 234 14]./255;
    ylim([min(result.vm.binnedXCoords{j}) max(result.vm.binnedXCoords{j})]);

    applyPaperFormatting
    %text(80,.85,'0.1 y units ','horiz','right','vert','middle');
    %%
    print([figureResultFolder filesep 'pdf' filesep 'Kip2Kymograph.pdf'], '-dpdf');
    export_fig([figureResultFolder filesep 'png' filesep 'Kip2Kymograph.png'], '-q101', '-r300');
