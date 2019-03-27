function [Xpk, Ypk] = findpeaks_stripped2(Yin)

pp1 = [Yin(1) Yin(1:end-1)];
pp2 = [Yin(2:end) Yin(end)];


Xpk = find((Yin > pp1) & (Yin > pp2));

if nargin > 1
    Ypk = Yin(Xpk);
end
