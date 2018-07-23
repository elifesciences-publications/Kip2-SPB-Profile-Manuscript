function [Ypk,Xpk] = findpeaks_stripped(Yin)
yIsRow = 1;
xIsRow = yIsRow;
y = Yin(:);
minH = -Inf;
refW = 'halfprom';
minT = 0;
minD = 0;
maxN = length(y);
x = (1:maxN)';

% find indices of all finite and infinite peaks and the inflection points
[iFinite,iInfinite,iInflect] = getAllPeaks(y);

% keep only the indices of finite peaks that meet the required
% minimum height and threshold
iPk = removePeaksBelowMinPeakHeight(y,iFinite,minH,refW);
iPk = removePeaksBelowThreshold(y,iPk,minT);


% find the indices of the largest peaks within the specified distance
idx = findPeaksSeparatedByMoreThanMinPeakDistance(y,x,iPk,minD);

% re-order and bound the number of peaks based upon the index vector
idx = keepAtMostNpPeaks(idx,maxN);

% use the index vector to fetch the correct peaks.
iPk = iPk(idx);

[Ypk,Xpk] = assignOutputs(y,x,iPk,yIsRow,xIsRow);
  

%--------------------------------------------------------------------------
function [iPk,iInf,iInflect] = getAllPeaks(y)
% fetch indices all infinite peaks
iInf = find(isinf(y) & y>0);

% temporarily remove all +Inf values
yTemp = y;
yTemp(iInf) = NaN;

% determine the peaks and inflection points of the signal
[iPk,iInflect] = findLocalMaxima(yTemp);


%--------------------------------------------------------------------------
function [iPk, iInflect] = findLocalMaxima(yTemp)
% bookend Y by NaN and make index vector
yTemp = [NaN; yTemp; NaN];
iTemp = (1:length(yTemp)).';

% keep only the first of any adjacent pairs of equal values (including NaN).
yFinite = ~isnan(yTemp);
iNeq = [1; 1 + find((yTemp(1:end-1) ~= yTemp(2:end)) & ...
                    (yFinite(1:end-1) | yFinite(2:end)))];
iTemp = iTemp(iNeq);

% take the sign of the first sample derivative
s = sign(diff(yTemp(iTemp)));

% find local maxima
iMax = 1 + find(diff(s)<0);

% find all transitions from rising to falling or to NaN
iAny = 1 + find(s(1:end-1)~=s(2:end));

% index into the original index vector without the NaN bookend.
iInflect = iTemp(iAny)-1;
iPk = iTemp(iMax)-1;


%--------------------------------------------------------------------------
function iPk = removePeaksBelowMinPeakHeight(Y,iPk,Ph,widthRef)
if ~isempty(iPk) 
  iPk = iPk(Y(iPk) > Ph);
  if isempty(iPk) && ~strcmp(widthRef,'halfheight')
    warning(message('signal:findpeaks:largeMinPeakHeight', 'MinPeakHeight', 'MinPeakHeight'));
  end
end
    
%--------------------------------------------------------------------------
function iPk = removePeaksBelowThreshold(Y,iPk,Th)
base = max(Y(iPk-1),Y(iPk+1));
iPk = iPk(Y(iPk)-base >= Th);


%--------------------------------------------------------------------------
function [iPkOut,bPk,bxPk,byPk,wxPk] = combinePeaks(iPk,iInf)
iPkOut = union(iPk,iInf);
bPk = zeros(0,1);
bxPk = zeros(0,2);
byPk = zeros(0,2);
wxPk = zeros(0,2);

%--------------------------------------------------------------------------
function idx = findPeaksSeparatedByMoreThanMinPeakDistance(y,x,iPk,Pd)
% Start with the larger peaks to make sure we don't accidentally keep a
% small peak and remove a large peak in its neighborhood. 

if isempty(iPk) || Pd==0
  IONE = ones('like',getIZERO);
  idx = (IONE:length(iPk)).';
  return
end

% copy peak values and locations to a temporary place
pks = y(iPk);
locs = x(iPk);

% Order peaks from large to small
if coder.target('MATLAB')
    [~, sortIdx] = sort(pks,'descend');
else
    sortIdx = coder.internal.sortIdx(pks,'d');
end
locs_temp = locs(sortIdx);

idelete = ones(size(locs_temp))<0;
for i = 1:length(locs_temp)
  if ~idelete(i)
    % If the peak is not in the neighborhood of a larger peak, find
    % secondary peaks to eliminate.
    idelete = idelete | (locs_temp>=locs_temp(i)-Pd)&(locs_temp<=locs_temp(i)+Pd); 
    idelete(i) = 0; % Keep current peak
  end
end

% report back indices in consecutive order
idx = sort(sortIdx(~idelete));





%--------------------------------------------------------------------------
function idx = keepAtMostNpPeaks(idx,Np)

if length(idx)>Np
  idx = idx(1:Np);
end


%--------------------------------------------------------------------------
function [YpkOut,XpkOut] = assignOutputs(y,x,iPk,yIsRow,xIsRow)
coder.internal.prefer_const(yIsRow,xIsRow);

% fetch the coordinates of the peak
Ypk = y(iPk);
Xpk = x(iPk);

% preserve orientation of Y
if yIsRow
  YpkOut = Ypk.';
else
  YpkOut = Ypk;
end

% preserve orientation of X
if xIsRow
  XpkOut = Xpk.';
else
  XpkOut = Xpk;
end


%--------------------------------------------------------------------------
function y = getIZERO
% Return zero of the indexing type: double 0 in MATLAB,
% coder.internal.indexInt(0) for code generation targets.
if coder.target('MATLAB')
    y = 0;
else
    y = coder.internal.indexInt(0);
end

% [EOF]
