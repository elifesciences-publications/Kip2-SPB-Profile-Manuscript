function compareConditions()
% compareConditions: Generates mean / 95% CI plots of different experimental conditions / bins
% used in the manuscript
%
%   Usage:
%   compareConditions()

% Â© 2018, ETH Zurich, Lukas Widmer (l.widmer@gmail.com)

load(['..' filesep 'analyzedData' filesep 'binnedProfiles.mat'], 'conditionResults');

%%
[conditions, maxLength] = displayConditions(conditionResults);

%%
conditionCombinations = {};
% 
% % wt bud data comparison
% for currentConditionIndex = 0:6
%     combo = struct;
%     combo.plotRegressionLine = false;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = currentConditionIndex + fliplr([43 62 326]);
%     %combo.adjustFluorescenceFromCondition = currentConditionIndex + [43 62 326];
%     %combo.adjustFluorescenceToReferenceCondition = currentConditionIndex + [43 43 43];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (bud), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 1.5]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.adjustFluorescenceFromCondition = [1 3 35];
%     combo.adjustFluorescenceToReferenceCondition = [1 1 1];
%     conditionCombinations{end+1} = combo;
% end
% 
% % wt mom data comparison
% for currentConditionIndex = 0:6
%     combo = struct;
%     combo.plotRegressionLine = false;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = currentConditionIndex + fliplr([53 71 336]);
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (mom), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 1.5]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.adjustFluorescenceFromCondition = [1 3 35];
%     combo.adjustFluorescenceToReferenceCondition = [1 1 1];
%     conditionCombinations{end+1} = combo;
% end


% 15100 vs 15100b, bud
for currentConditionIndex = [43 62 326]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
    %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    combo.slopeTag = ['wt' int2str(currentConditionIndex)];
    conditionCombinations{end+1} = combo;
end

i = 1;
adustVector = [1 3 35];
for currentConditionIndex = [43 62 326]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* adustVector(i);
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    i = i + 1;
end


% 15100 vs 15100b, mom
for currentConditionIndex = [53 71 336]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
    %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (mom), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end

i = 1;
adustVector = [1 3 35];
for currentConditionIndex = [53 71 336]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* adustVector(i);
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (mom), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    i = i + 1;
end

% bfa1Del bud
for currentConditionIndex = [344]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta (bud), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    combo.adjustFluorescenceFromCondition = 35*[1 1 1 1 1 1];
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];    
    conditionCombinations{end+1} = combo;
    
    combo.restrictToDistalIndex = -38;
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta (bud, no mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;     
    
    combo.restrictToDistalIndex = 38;
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta (bud, mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;   
end

% bub2Del bud
for currentConditionIndex = [363]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP bub2\Delta (bud), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    combo.adjustFluorescenceFromCondition = 35*[1 1 1 1 1 1];
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];    
    conditionCombinations{end+1} = combo;
    
    combo.restrictToDistalIndex = -40;
    combo.titleText = 'Kip2-3xsfGFP bub2\Delta (bud, no mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;     
    
    combo.restrictToDistalIndex = 40;
    combo.titleText = 'Kip2-3xsfGFP bub2\Delta (bud, mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;   
end

% bfa1bub2Del bud
for currentConditionIndex = [383]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta bub2\Delta (bud), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    combo.adjustFluorescenceFromCondition = 35*[1 1 1 1 1 1];
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];   
    combo.plotSlopeTag = 'wt43';
    conditionCombinations{end+1} = combo;
    
    combo.restrictToDistalIndex = -42;
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta bub2\Delta (bud, no mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;     
    
    combo.restrictToDistalIndex = 42;
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta bub2\Delta (bud, mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;   
end


% bfa1Del mom
for currentConditionIndex = [353]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta (mom), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    combo.adjustFluorescenceFromCondition = 35*[1 1 1 1 1 1];
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];    
    conditionCombinations{end+1} = combo;
end

% bub2Del mom
for currentConditionIndex = [373]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP bub2\Delta (mom), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    combo.adjustFluorescenceFromCondition = 35*[1 1 1 1 1 1];
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];    
    conditionCombinations{end+1} = combo;
end

% bfa1bub2Del mom
for currentConditionIndex = [392]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP bfa1\Delta bub2\Delta (mom), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
    combo.adjustFluorescenceFromCondition = 35*[1 1 1 1 1 1];
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];  
    conditionCombinations{end+1} = combo;
end


%%% TODO

% wt different strains / different days
% for currentConditionIndex = 0:6
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = currentConditionIndex + [102 120 54 178 195];
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (bud), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% 
% end
% 
% for currentConditionIndex = 0:5
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = currentConditionIndex + [45 63 187 205];
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (mom), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end

% 
% % 15100 vs 15100b, bud
% for currentConditionIndex = [178 195]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (bud), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% % 15100 vs 15100b, mom
% for currentConditionIndex = [187 205]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (mom), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% % 15692 vs 15692b (bud)
% for currentConditionIndex = [220 239]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP bfa1\delta (bud), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% % 15692 vs 15692b (mom)
% for currentConditionIndex = [229 212]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP bfa1\delta (mom), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% % 15693 vs 15693b (bud)
% for currentConditionIndex = [248 267]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP bub2\delta (bud), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% % 15693 vs 15693b (mom)
% for currentConditionIndex = [258 276]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP bub2\delta (mom), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% % 15794 vs 15794b (bud)
% for currentConditionIndex = [286 305]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP bfa1\Delta bub2\Delta (bud), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% % 15794 vs 15794b (mom)
% for currentConditionIndex = [295]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     %combo.adjustFluorescenceFromCondition = [9 11 3 19 21];
%     %combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1];
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP bfa1\Delta bub2\Delta (mom), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% end
% for currentConditionIndex = [9]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; 
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (bud), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% 
%     combo.restrictToDistalIndex = -2;
%     combo.titleText = 'Kip2-3xsfGFP (bud, no mom mt), mean profile / 95% CI';
%     combo.legendLocation = 'nw'; 
%     combo.legendBuffer = [5 -5];
%     combo.ylimsGreen = [-0.05 0.9]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.restrictToDistalIndex = 2;
%     combo.titleText = 'Kip2-3xsfGFP (bud, mom mt), mean profile / 95% CI';
%     combo.legendLocation = 'nw'; 
%     combo.legendBuffer = [5 -5];
%     combo.ylimsGreen = [-0.05 0.9]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.alignAt = 'both';
%     combo.titleText = 'Kip2-3xsfGFP (bud), mean profile / 95% CI';
%     combo.restrictToDistalIndex = 0;
%     combo.legendLocation = 'ne'; 
%     combo.legendBuffer = [-10 -10];
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% 
%     combo.restrictToDistalIndex = -2;
%     combo.titleText = 'Kip2-3xsfGFP (bud, no mom mt), mean profile / 95% CI';
%     combo.legendLocation = 'nw'; 
%     combo.legendBuffer = [5 -5];
%     combo.ylimsGreen = [-0.05 0.9]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.restrictToDistalIndex = 2;
%     combo.titleText = 'Kip2-3xsfGFP (bud, mom mt), mean profile / 95% CI';
%     combo.legendLocation = 'nw'; 
%     combo.legendBuffer = [5 -5];
%     combo.ylimsGreen = [-0.05 0.9]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
% %NEW Kip2 wt metaphase (distal)
% for currentConditionIndex = [19]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+4):-1:currentConditionIndex;
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     combo.legendLocation = 'n'; %NorthEastOutside
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip2-3xsfGFP (distal), mean profile / 95% CI';
%     combo.grayoutSpindleArea = false;
%     combo.alignAt = 'bothSPB';
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.legendLocation = 'nw'; %NorthEastOutside
%     combo.legendBuffer = [5 -5];    
%     combo.ylimsGreen = [-0.05 0.9]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.legendLocation = 'ne'; %NorthEastOutside
%     combo.legendBuffer = [-10 -10];
%     combo.alignAt = 'both';
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
% 
%     combo.legendLocation = 'nw'; %NorthEastOutside
%     combo.legendBuffer = [5 -5];
%     combo.ylimsGreen = [-0.05 0.9]*10;
%     conditionCombinations{end+1} = combo;
% end
% 
%  
%NEW S63A
for currentConditionIndex = [79]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+8):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'nw';
    combo.legendBuffer = [10 -10];
    combo.titleText = 'Kip2-S63A-3xsfGFP (bud), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.xlims    = [-3.7 0.8];
    combo.xlimsSPB = [-0.5 4];
    combo.ylimsGreen = [-0.1 1.5]*10;
    combo.alignAt = 'bothSPB';
    conditionCombinations{end+1} = combo;
    combo.xAspect = 1.25; %1.25
    conditionCombinations{end+1} = combo;
    
    combo.xAspect = 1; %1.25
    combo.adjustFluorescenceFromCondition = 3 * ones(1, length(combo.conditions));
    combo.adjustFluorescenceToReferenceCondition = ones(1, length(combo.conditions));
    combo.plotSlopeTag = 'wt43';
    conditionCombinations{end+1} = combo;
    combo.xAspect = 1.25; %1.25
    conditionCombinations{end+1} = combo;
    
    combo.xAspect = 1;
    combo.restrictToDistalIndex = -6;
    combo.titleText = 'Kip2-S63A-3xsfGFP (bud, no mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;     
    
    combo.restrictToDistalIndex = 6;
    combo.titleText = 'Kip2-S63A-3xsfGFP (bud, mom mt), mean profile / 95% CI';
    conditionCombinations{end+1} = combo;         
    
    %combo.restrictToDistalIndex = 0;
%     combo.titleText = 'Kip2-S63A-3xsfGFP (bud), mean profile / 95% CI';
%     combo.xAspect = 1;
%     combo.alignAt = 'both';
%     conditionCombinations{end+1} = combo;
%     combo.xAspect = 1.25;
%     conditionCombinations{end+1} = combo;
end

%NEW S63A distal
for currentConditionIndex = [89]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+8):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    %combo.alignAt = 'plusEnd';
    combo.legendLocation = 'nw';
    combo.legendBuffer = [10 -10];
    combo.titleText = 'Kip2-S63A-3xsfGFP (distal), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.xlims    = [-3.7 0.8];
    combo.xlimsSPB = [-0.5 4];
    combo.ylimsGreen = [-0.1 1.5]*10;
    %conditionCombinations{end+1} = combo;
    combo.alignAt = 'bothSPB';
    conditionCombinations{end+1} = combo;
    combo.xAspect = 1.25;
    conditionCombinations{end+1} = combo;
    
    combo.xAspect = 1; %1.25
    combo.adjustFluorescenceFromCondition = 3 * ones(1, length(combo.conditions));
    combo.adjustFluorescenceToReferenceCondition = ones(1, length(combo.conditions));
    conditionCombinations{end+1} = combo;
    combo.xAspect = 1.25; %1.25
    conditionCombinations{end+1} = combo;    
%     combo.alignAt = 'both';
%     conditionCombinations{end+1} = combo;
%     combo.xAspect = 1.25;
%     conditionCombinations{end+1} = combo;
end
% 

%NEW Kip3
for currentConditionIndex = [98]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    %combo.alignAt = 'plusEnd';
    combo.legendLocation = 'n';
    combo.legendBuffer = [-50 -10];
    combo.titleText = 'Kip3-3xsfGFP (bud), mean profile / 95% CI';
    combo.grayoutSpindleArea = true;
    %conditionCombinations{end+1} = combo;
    combo.alignAt = 'bothSPB';
    combo.ylimsGreen = [-0.1 2]*10;
    combo.yTicks = [0 5 10 15 20];
    conditionCombinations{end+1} = combo;
    
%     combo.legendBuffer = [10 -10];
%     combo.alignAt = 'both';
%     conditionCombinations{end+1} = combo;
end

%NEW Kip3 distal
for currentConditionIndex = [106]
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 0;
    combo.conditions = (currentConditionIndex+1):-1:currentConditionIndex;
    combo.legendConditions = fliplr(combo.conditions);
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    %combo.alignAt = 'plusEnd';
    combo.legendLocation = 'n';
    combo.legendBuffer = [-50 -10];
    combo.titleText = 'Kip3-3xsfGFP (mom), mean profile / 95% CI';
    combo.grayoutSpindleArea = true;
    %conditionCombinations{end+1} = combo;
    combo.alignAt = 'bothSPB';
    combo.ylimsGreen = [-0.1 2]*10;
    combo.yTicks = [0 5 10 15 20];
    conditionCombinations{end+1} = combo;
    
%     combo.legendBuffer = [10 -10];
%     combo.alignAt = 'both';
%     conditionCombinations{end+1} = combo;
end

% wt Jan, bud (no mom MT)
for currentConditionIndex = 110
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = -10;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* 1;
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud, no mom MT), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end

% wt Feb, bud (no mom MT)
for currentConditionIndex = 62
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = -4;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* 3;
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud, no mom MT), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end


% wt Dec, bud (no mom MT)
for currentConditionIndex = 326
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = -36;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* 35;
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud, no mom MT), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end


% wt Jan, bud (mom MT)
for currentConditionIndex = 110
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 10;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* 1;
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud, mom MT), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end

% wt Feb, bud (mom MT)
for currentConditionIndex = 62
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 4;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* 3;
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud, mom MT), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end

% wt Dec, bud (mom MT)
for currentConditionIndex = 326
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 36;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* 35;
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (bud, mom MT), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end

% wt Jan, mom
for currentConditionIndex = 119
    combo = struct;
    combo.plotRegressionLine = true;
    combo.restrictToDistalIndex = 10;
    combo.conditions = (currentConditionIndex+5):-1:currentConditionIndex;
    combo.adjustFluorescenceFromCondition = [1 1 1 1 1 1] .* 1;
    combo.adjustFluorescenceToReferenceCondition = [1 1 1 1 1 1];
    combo.legendConditions = fliplr(combo.conditions);
    combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
    combo.legendLocation = 'n'; 
    combo.legendBuffer = [0 -10];
    combo.titleText = 'Kip2-3xsfGFP (mom), mean profile / 95% CI';
    combo.grayoutSpindleArea = false;
    combo.alignAt = 'bothSPB';
    combo.xlims    = [-2.4 0.6];
    combo.xlimsSPB = combo.xlims + 2.1;
    combo.ylimsGreen = [-0.1 1.5]*10;
    conditionCombinations{end+1} = combo;
end

% 
% %NEW Kip3 distal
% for currentConditionIndex = [72]
%     combo = struct;
%     combo.restrictToDistalIndex = 0;
%     combo.conditions = (currentConditionIndex+1):-1:currentConditionIndex;
%     combo.legendConditions = fliplr(combo.conditions);
%     combo.xlims    = [-2.4 0.6];
%     combo.xlimsSPB = combo.xlims + 2.1;
%     combo.yAxisLabel = 'GFP fluorescence along aMT (a.u.)';
%     %combo.alignAt = 'plusEnd';
%     combo.legendLocation = 'n';
%     combo.legendBuffer = [0 -10];
%     combo.titleText = 'Kip3-3xsfGFP (distal), mean profile / 95% CI';
%     combo.grayoutSpindleArea = true;
%     %conditionCombinations{end+1} = combo;
%     combo.alignAt = 'bothSPB';
%     combo.ylimsGreen = [-0.1 2]*10;
%     conditionCombinations{end+1} = combo;
%     
%     combo.legendBuffer = [10 -10];
%     combo.legendLocation = 'nw';
%     combo.alignAt = 'both';
%     conditionCombinations{end+1} = combo;
% end


%%
figureFolder    = ['..' filesep 'figures' filesep 'profileComparisons'];
figureFolderSep = [figureFolder filesep];

if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

if ~exist([figureFolderSep 'pdf'], 'dir')
    mkdir([figureFolderSep 'pdf']);
end

if ~exist([figureFolderSep 'png'], 'dir')
    mkdir([figureFolderSep 'png']);
end

exportPlots = true;
applyBackgroundCorrection = 'bulk';
plotRed = false;
plotData = false;
debugFluoMapping = false;

if strcmp(applyBackgroundCorrection, 'none')
    for combinationIndex = 1:length(conditionCombinations)
        conditionCombinations{combinationIndex}.ylimsGreen = conditionCombinations{combinationIndex}.ylimsGreen + [1.1 2]*10;
    end
end
ylimsRed = [0.2 1.5]*1e5;

minimumCIprofiles = 8;
lineWidth = 1.181;
opacity = 0.1;
showLegend = true;

plotError = 'CI';

fileSuffix = '';

fileSuffix = [fileSuffix '-bgCorrected-' applyBackgroundCorrection '-' plotError];

if plotRed
    fileSuffix = [fileSuffix '-red'];
end

maxIndex = round(maxLength ./ (4/30));

savedSlopes = struct();
linearModelCache = struct();

%%
for combinationIndex = 1:length(conditionCombinations)
    combo = conditionCombinations{combinationIndex};
    conditionsToPlot = combo.conditions;
    plusEndLocations  = nan(1, length(combo.legendConditions));
    meanPlusEndValues = nan(1, length(combo.legendConditions));
    semPlusEndValues  = nan(1, length(combo.legendConditions));
    legendText = conditions(combo.legendConditions);
    for currentConditionIndex = 1:length(combo.legendConditions)
        c = conditionResults{combo.legendConditions(currentConditionIndex)};
        legendText{currentConditionIndex} = regexprep(legendText{currentConditionIndex}, '^(.+bin-)(.*?)-(.*)$', '$2 - $3');
        if combo.restrictToDistalIndex ~= 0
            if contains(c.condition, 'Distal')
                % Distal MTs always have proximal ones due to our analysis
            else
                if combo.restrictToDistalIndex < 0
                    rIndex = -combo.restrictToDistalIndex;
                    restrictedIndices = find(~contains(c.cellNames, conditionResults{rIndex}.cellNames));
                else
                    restrictedIndices = find(contains(c.cellNames, conditionResults{combo.restrictToDistalIndex}.cellNames));
                end
                c.greenIntensities = c.greenIntensities(:, restrictedIndices);
            end
        end
        legendText{currentConditionIndex} = [legendText{currentConditionIndex} sprintf(' (n = %3i)', size(c.greenIntensities,2))];


    end
    currentFigure = figure();
    hold on
    colorOrder = get(gca, 'ColorOrder');
    nCond = 0;
    handles = [];
    
    if combo.grayoutSpindleArea
        switch combo.alignAt
            case 'bothSPB'
                if combo.grayoutSpindleArea
                    %grayArea = area([0 max(xlimsSPB)], max(combo.ylimsGreen) * [1 1], min(combo.ylimsGreen), 'FaceColor', [1 1 1]*0.5, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                    grayArea = patch([0 min(combo.xlimsSPB) min(combo.xlimsSPB) 0], [min(combo.ylimsGreen) * [1 1] max(combo.ylimsGreen) * [1 1]], [1 1 1]*0.5, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
            case 'SPB'
                if combo.grayoutSpindleArea
                    grayArea = area([0 max(combo.xlimsSPB)], max(combo.ylimsGreen) * [1 1], min(combo.ylimsGreen), 'FaceColor', [1 1 1]*0.5, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
                end
        end
    end
    
    for currentConditionIndex = conditionsToPlot
        nCond = find(combo.legendConditions == currentConditionIndex);

        if nCond <= 7
            style = '';
        else
            style = '--';
        end

        c = conditionResults{currentConditionIndex};

        nMicrotubules = size(c.greenIntensities, 2);
        
        switch combo.alignAt
            case 'plusEnd'
                greenIntensitiesShifted = c.greenIntensitiesShifted;
                redIntensitiesShifted   = c.redIntensitiesShifted;
                shiftedLocations        = c.lengthVectorShifted;
                currentXlims            = combo.xlims;
                currentXLabel           = 'Distance from plus end ({\mu}m)';
                SPBlocation             = mean(c.lengthsUM);
            case 'SPB'
                greenIntensitiesShifted = c.greenIntensitiesSPBShifted;
                redIntensitiesShifted   = c.redIntensitiesSPBShifted;
                shiftedLocations        = c.lengthVectorSPBShifted;
                currentXlims            = combo.xlimsSPB;
                currentXLabel           = 'Distance from SPB ({\mu}m)';
                SPBlocation             = 0;
            case 'both'
                greenIntensitiesShifted = c.greenIntensitiesCenterShifted;
                redIntensitiesShifted   = c.redIntensitiesCenterShifted;
                shiftedLocations        = -c.lengthVectorCenterShifted;
                currentXlims            = combo.xlims;
                currentXLabel           = 'Distance from plus end ({\mu}m)';
                SPBlocation             = (c.meanMeanOffset-0.5)*(4/30);
            case 'bothSPB'
                greenIntensitiesShifted = c.greenIntensitiesCenterShifted;
                redIntensitiesShifted   = c.redIntensitiesCenterShifted;
                shiftedLocations        = c.lengthVectorCenterShifted;
                currentXlims            = combo.xlimsSPB;
                currentXLabel           = 'Distance from SPB ({\mu}m)';
                SPBlocation             = (c.meanMeanOffset-0.5)*(4/30);
                shiftedLocations        = -(shiftedLocations - SPBlocation);
                plusEndLocation         = SPBlocation + 2/30;
                SPBlocation             = 0;
                
            case default
                error('Unknown alignment option');
        end

        if isfield(combo, 'adjustFluorescenceToReferenceCondition')
            referenceCondition = combo.adjustFluorescenceToReferenceCondition(nCond);
            fromCondition = combo.adjustFluorescenceFromCondition(nCond);
            %currentConditionIndex
            
            cFrom = conditionResults{fromCondition};
            cRef = conditionResults{referenceCondition};
            
            cacheTag = ['from' int2str(fromCondition) 'to' int2str(referenceCondition)];
            
            if isfield(linearModelCache, cacheTag)
            	fitresult = linearModelCache.(cacheTag).fitresult;
                gof = linearModelCache.(cacheTag).gof;
            else
                fprintf('Current   condition: %i, \t%s\n', currentConditionIndex, c.condition);
                fprintf('From      condition: %i, \t%s\n', fromCondition, cFrom.condition);
                fprintf('Reference condition: %i, \t%s\n', referenceCondition, cRef.condition);

                quantilesForCalibration = [0.005:0.002:0.995];

                nMicrotubulesRef = size(cRef.greenIntensities, 2);
                refIntensities = [];

                refPlusEndLoc = cRef.offset;

                for j = 1:nMicrotubulesRef
                    if cRef.lengths(j) ~= 0
                        refSPBloc = refPlusEndLoc + cRef.lengths(j);
                        currentMtIntensities = cRef.greenIntensitiesShifted(refPlusEndLoc:refSPBloc, j);
                        refIntensities = [refIntensities; currentMtIntensities];
                    end
                end

                assert(all(isfinite(refIntensities)));

                nMicrotubulesFrom = size(cFrom.greenIntensities, 2);
                fromIntensities = [];

                fromPlusEndLoc = cFrom.offset;

                for j = 1:nMicrotubulesFrom
                    if cFrom.lengths(j) ~= 0
                        fromSPBloc = fromPlusEndLoc + cFrom.lengths(j);
                        currentMtIntensities = cFrom.greenIntensitiesShifted(fromPlusEndLoc:fromSPBloc, j);
                        fromIntensities = [fromIntensities; currentMtIntensities];
                    end
                end            

                assert(all(isfinite(fromIntensities)));

                allQuantiles = quantile(fromIntensities, quantilesForCalibration);
                allQuantilesRef = quantile(refIntensities, quantilesForCalibration);

                [fitresult, gof] = powerModelFit(allQuantiles, allQuantilesRef, debugFluoMapping)

                linearModelCache.(cacheTag).fitresult = fitresult;
                linearModelCache.(cacheTag).gof = gof;
            end
            
            originalSize = size(greenIntensitiesShifted);
            greenIntensitiesShifted = fitresult(greenIntensitiesShifted);
            greenIntensitiesShifted = reshape(greenIntensitiesShifted, originalSize);

            originalSize = size(c.backgroundIntensities);
            c.backgroundIntensities = fitresult(c.backgroundIntensities);
            c.backgroundIntensities = reshape(c.backgroundIntensities, originalSize);            
            if debugFluoMapping
                % Debug plot
                
                figure;
                hold on;
                %bgHandle = plot(bgQuantiles, bgQuantilesRef, '.--');
                %spbHandle = plot(spbQuantiles, spbQuantilesRef, '.--');
                %plusEndHandle = plot(plusEndQuantiles, plusEndQuantilesRef, '.--');
                allHandle = plot(allQuantiles, allQuantilesRef, '.--');
                xlabel('Fluorescence (AU)');
                ylabel('Reference fluorescence (AU)');
                title(['Condition ' int2str(fromCondition) ' VS reference condition ' int2str(referenceCondition)]);
                %legend([bgHandle, spbHandle, plusEndHandle, allHandle], {'Background', 'SPB', 'Plus end', 'All'});
                legend([allHandle], {'Profile fluorescence Q-Q'});
                axis equal;
                
                
                
                figure(currentFigure);

            end
                
            
        end


        switch applyBackgroundCorrection
            case 'bulk'
                greenIntensitiesShifted = greenIntensitiesShifted - mean(c.backgroundIntensities);
            case 'bulkMedian'
                greenIntensitiesShifted = greenIntensitiesShifted - median(c.backgroundIntensities);
            case 'individual'
                
                if strcmp(c.dataSet.backgroundComputation, 'both')
                    bgIntensities = min(c.backgroundIntensities(1:2:end), c.backgroundIntensities(2:2:end));
                else
                    bgIntensities = c.backgroundIntensities;
                end
                greenIntensitiesShifted = greenIntensitiesShifted - repmat(bgIntensities',size(greenIntensitiesShifted,1),1);
            case 'none'
                
            case default
                error('Unknown background correction type');
        end
     
        if combo.restrictToDistalIndex ~= 0
            if contains(c.condition, 'Distal')
                % Distal MTs always have proximal ones due to our analysis
            else
                if combo.restrictToDistalIndex < 0
                    rIndex = -combo.restrictToDistalIndex;
                    restrictedIndices = find(~contains(c.cellNames, conditionResults{rIndex}.cellNames));
                else
                    restrictedIndices = find(contains(c.cellNames, conditionResults{combo.restrictToDistalIndex}.cellNames));
                end
                greenIntensitiesShifted = greenIntensitiesShifted(:, restrictedIndices);
                c.lengthsUM = c.lengthsUM(restrictedIndices);
                nMicrotubules = length(restrictedIndices);

                %fprintf('N = %i\n', length(restrictedIndices))
            end
        end
        greenIntensitiesShifted = greenIntensitiesShifted * 1e-3; % Re-scale arbitrary units
        meanGreenIntensityCenterShifted = mean(greenIntensitiesShifted,2, 'omitnan');


        stdGreenIntensityCenterShifted = std(greenIntensitiesShifted,0,2, 'omitnan');
        nMTsCenterShifted = sum(~isnan(greenIntensitiesShifted),2);
        semGreenIntensityCenterShifted = stdGreenIntensityCenterShifted./sqrt(nMTsCenterShifted);

        semGreenIntensityCenterShifted(nMTsCenterShifted < minimumCIprofiles) = NaN;
        meanGreenIntensityCenterShifted(nMTsCenterShifted < minimumCIprofiles) = NaN;

        noNan = ~any(isnan(greenIntensitiesShifted),2)';

        if strcmp(plotError, 'CI')    
            meanCIfunc  = @(i) mean(i, 'omitnan');

            greenMeanCI = nan(2, length(shiftedLocations));
            availableProfiles = ~isnan(greenIntensitiesShifted);
            nAvailableProfiles = sum(availableProfiles, 2);
            indices = find(nAvailableProfiles >= minimumCIprofiles)';
            parfor j = indices
                currentIntensities = greenIntensitiesShifted(j, :);
                currentIntensities = currentIntensities(isfinite(currentIntensities));
                greenMeanCI(:, j) = bootci(5000, meanCIfunc, currentIntensities);
            end

            greenMeanCI(1, :) = meanGreenIntensityCenterShifted' - greenMeanCI(1, :);
            greenMeanCI(2, :) = greenMeanCI(2, :) - meanGreenIntensityCenterShifted';
        end


        if plotRed
            meanRedIntensityCenterShifted = mean(redIntensitiesShifted, 2, 'omitnan');
            stdRedIntensityCenterShifted = std(redIntensitiesShifted, 0, 2, 'omitnan');
            semRedIntensityCenterShifted = stdRedIntensityCenterShifted./sqrt(nMTsCenterShifted);

            semRedIntensityCenterShifted(nMTsCenterShifted < minimumCIprofiles) = NaN;
            meanRedIntensityCenterShifted(nMTsCenterShifted < minimumCIprofiles) = NaN;

            %redMeanCInoNan = bootci(10000, meanCIfunc, c.redIntensitiesCenterShifted(noNan, :)');

            if strcmp(plotError, 'CI')                
                redMeanCI = nan(2, length(shiftedLocations));
                availableProfiles = ~isnan(redIntensitiesShifted);
                nAvailableProfiles = sum(availableProfiles, 2);
                indices = find(nAvailableProfiles >= 10)';
                parfor j = indices
                    currentIntensities = redIntensitiesShifted(j, :);
                    currentIntensities = currentIntensities(isfinite(currentIntensities));
                    redMeanCI(:, j) = bootci(5000, meanCIfunc, currentIntensities);
                end

                redMeanCI(1, :) = meanRedIntensityCenterShifted' - redMeanCI(1, :);
                redMeanCI(2, :) = redMeanCI(2, :) - meanRedIntensityCenterShifted';
            end
        end



        currentColor = colorOrder(mod(nCond-1,7)+1,:);

        if plotData
            for j = 1:nMicrotubules
                if c.lengths(j) == 0
                    if plotRed
                        subplot(2,1,1);
                    end
                    hold on;
                    plot(shiftedLocations, greenIntensitiesShifted(:, j), 'y');
                else
                    if plotRed
                        subplot(2,1,1);
                    end
                    hold on;
                    foo = plot(shiftedLocations, greenIntensitiesShifted(:, j), style, 'Color', currentColor);
                    foo.Color(4) = opacity;

                    if plotRed
                        subplot(2,1,2);
                        hold on;
                        foo = plot(shiftedLocations, redIntensitiesShifted(:, j), style, 'Color', currentColor);
                        foo.Color(4) = opacity;
                    end
                end
            end
        end
        if plotRed
            subplot(2,1,1);
        end

        switch plotError
            case 'CI'
                H = shadedErrorBar(shiftedLocations, meanGreenIntensityCenterShifted, greenMeanCI, 'lineProps', {style, 'Color', currentColor, 'LineWidth', lineWidth});
            case 'SEM'
                H = shadedErrorBar(shiftedLocations, meanGreenIntensityCenterShifted, semGreenIntensityCenterShifted, 'lineProps', {style, 'Color', currentColor, 'LineWidth', lineWidth});
            case 'STD'
                H = shadedErrorBar(shiftedLocations, meanGreenIntensityCenterShifted, stdGreenIntensityCenterShifted, 'lineProps', {style, 'Color', currentColor, 'LineWidth', lineWidth});
            case ''
                % No error bar
            case default
                
        end
        %H = shadedErrorBar(c.lengthVectorCenterShifted, meanGreenIntensityCenterShifted, semGreenIntensityCenterShifted, 'lineProps', {style, 'Color', currentColor, 'LineWidth', 2});
        hold on
        switch combo.alignAt
            case 'both'
                H_SPB = plot(-[SPBlocation SPBlocation], combo.ylimsGreen,  '--', 'Color', currentColor);
            case 'plusEnd'
                H_SPB = plot(-[SPBlocation SPBlocation], combo.ylimsGreen,  '--', 'Color', currentColor);
        end
        
        if plotRed
            subplot(2,1,2);
            
            switch plotError
                case 'CI'
                    H = shadedErrorBar(shiftedLocations, meanRedIntensityCenterShifted, redMeanCI, 'lineProps', {style, 'Color', currentColor, 'LineWidth', lineWidth});
                case 'SEM'
                    H = shadedErrorBar(shiftedLocations, meanRedIntensityCenterShifted, semRedIntensityCenterShifted, 'lineProps', {style, 'Color', currentColor, 'LineWidth', lineWidth});
                case 'STD'
                    H = shadedErrorBar(shiftedLocations, meanRedIntensityCenterShifted, stdRedIntensityCenterShifted, 'lineProps', {style, 'Color', currentColor, 'LineWidth', lineWidth});
                case ''
                    % No error bar
                case default
                    error('Unknown error type');
            end
    %        H = shadedErrorBar(c.lengthVectorCenterShifted, meanRedIntensityCenterShifted, semRedIntensityCenterShifted, 'lineProps', {style, 'Color', currentColor, 'LineWidth', 2});
            %H = shadedErrorBar(shiftedLocations, meanRedIntensityCenterShifted, redMeanCI, 'lineProps', {style, 'Color', currentColor, 'LineWidth', 2});

            hold on
        end
        
        [~, closestIndex] = min(abs(shiftedLocations - plusEndLocation));
        plusEndLocations(nCond)  = shiftedLocations(closestIndex);
        meanPlusEndValues(nCond) = meanGreenIntensityCenterShifted(closestIndex);
        semPlusEndValues(nCond)  = semGreenIntensityCenterShifted(closestIndex);
        handles(end+1) = H.mainLine;
        
    end
    %%
    if plotRed
        subplot(2,1,1);
    end
    if ~isempty(applyBackgroundCorrection)
        plot(currentXlims, [0 0], '--k');
    end
    

    
    xlim(currentXlims);
    ylim(combo.ylimsGreen);
    if isfield(combo, 'yTicks')
        set(gca, 'YTick', combo.yTicks);
    end
    %set(gca, 'XDir', 'reverse');
    ylabel(combo.yAxisLabel);
    

    title(combo.titleText);
    
    %set(gca, 'YScale', 'log');
    grid on;

    switch combo.alignAt
        case 'bothSPB'
            H_SPB = plot([0 0], combo.ylimsGreen,  '--', 'Color', [0 0 0]);
        case 'SPB'
            H_SPB = plot([0 0], combo.ylimsGreen,  '--', 'Color', [0 0 0]);
    end

    if plotRed
        subplot(2,1,2);
        xlim(currentXlims);
        ylim(ylimsRed);
        ylabel('Spc42-mCherry Fluorescence (AU)');
        grid on;
    end
    xlabel(currentXLabel);
%     [lgd,icons,plots,txt] = legend(fliplr(handles), legendText, 'Location', combo.legendLocation);
%     legendLines = icons((1+length(legendText)):2:end);
%     if length(legendLines) > 0
%          (legendLines(1).XData(2) - legendLines(1).XData(1))
%     end
%     for i=1:length(legendLines)
%         legendLines(i).XData(1) = legendLines(i).XData(1) - (legendLines(i).XData(2) - legendLines(i).XData(1)) * 0.5;
%     end
%     legend boxoff;
    
    if plotRed
        subplot(2,1,1);
    end
    foo2 = gca;
    foo2.XLabel.Color = [0 0 0];
    foo2.YLabel.Color = [0 0 0];
    foo2.XColor = [0 0 0];
    foo2.YColor = [0 0 0];
    
    if combo.plotRegressionLine

        linearModel = fitlm(plusEndLocations, meanPlusEndValues, 'Weights', (1./semPlusEndValues).^2);
        linearModelLineX = currentXlims';
        linearModelLineX(1) = 0;
        fitY = predict(linearModel, linearModelLineX);
        hFit = plot(linearModelLineX, fitY, '--', 'Color', [0.4 0.4 0.4 0.5], 'LineWidth', lineWidth);    
        %handles = [hFit handles];
        if plotRed
            yCoord = min(combo.ylimsGreen) - 0.3*mean(combo.ylimsGreen);
        else
            yCoord = 1;
        end
        text(mean(currentXlims), yCoord, sprintf('slope = %.3g ± %.3g, Y-intercept = %.3g ± %.3g', linearModel.Coefficients.Estimate(2), linearModel.Coefficients.SE(2), linearModel.Coefficients.Estimate(1), linearModel.Coefficients.SE(1)), 'HorizontalAlignment', 'center');
        fprintf('plot: %i, slope = %.3g ± %.3g, Y-intercept = %.3g ± %.3g\n', combinationIndex, linearModel.Coefficients.Estimate(2), linearModel.Coefficients.SE(2), linearModel.Coefficients.Estimate(1), linearModel.Coefficients.SE(1))
        %legendText{end+1} = sprintf('slope = %g ï¿½ %g, Y-intercept = %g ï¿½ %g', linearModel.Coefficients.Estimate(2), linearModel.Coefficients.SE(2), linearModel.Coefficients.Estimate(1), linearModel.Coefficients.SE(1))
        if isfield(combo, 'slopeTag')
            savedSlopes.(combo.slopeTag) = linearModel;
        end
        
        if isfield(combo, 'plotSlopeTag')
            fitY = predict(savedSlopes.(combo.plotSlopeTag), linearModelLineX);
            plot(linearModelLineX, fitY, '-.', 'Color', [0 0 0 0.3], 'LineWidth', lineWidth);   
        end
    end
    
    if plotRed
        subplot(2,1,2);
        foo2 = gca;
        foo2.XLabel.Color = [0 0 0];
        foo2.YLabel.Color = [0 0 0];
        foo2.XColor = [0 0 0];
        foo2.YColor = [0 0 0];
    end 
    
    applyPaperFormatting
    
    if showLegend
        [legend_h,object_h,plot_h,text_str] = legendflex(fliplr(handles), legendText, 'anchor', {combo.legendLocation, combo.legendLocation}, ...
        'buffer', combo.legendBuffer, ...
        'xscale', 0.4, ...
        'box', 'off', ...
        'title', 'Bin ({\mu}m)', ...
        'FontSize', 10);
        legend_h.Color = [1 1 1 0.8];
    end
    

    
    foo = gcf;
    foo.Color = 'white';
    if plotRed
        foo.Position(4) = foo.Position(4)*2.42;
    end
    if isfield(combo, 'xAspect')
        foo.Position(3) = foo.Position(3)*combo.xAspect;
    end
    if contains(combo.legendLocation, 'Outside', 'IgnoreCase', true)
        foo.Position(3) = foo.Position(3)*1.4;
    end
    
   
    
    movegui(foo, 'center')
    %%
    currentFileSuffix = fileSuffix;
    currentFileSuffix = [currentFileSuffix '-' combo.alignAt];
    if isfield(combo, 'adjustFluorescenceToReferenceCondition')
        currentFileSuffix = [currentFileSuffix '-refNormalized'];
    end
    
    if exist('grayArea', 'var') && any(grayArea == get(gca,'children'))
        uistack(grayArea,'bottom')
    end
    drawnow
    if exportPlots
        print([figureFolderSep 'pdf' filesep 'paper-' int2str(combinationIndex) '-cond' char(strjoin(string(conditionsToPlot), '_')) currentFileSuffix '-print-dpdf.pdf'], '-dpdf');
        export_fig([figureFolderSep 'png' filesep 'paper-' int2str(combinationIndex) '-cond' char(strjoin(string(conditionsToPlot), '_')) currentFileSuffix '.png'], '-q101', '-r300');
    else
        pause
    end
    
    close(foo);
end