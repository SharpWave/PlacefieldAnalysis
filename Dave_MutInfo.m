function [] = Dave_MutInfo(PosVec,NumPosBins,ActVec,NumActBins)
% Dave_MutInfo(PosVec,ActVec,GoodPoints)
% PosVec is a list of positions (bin # as a function of time)
% ActVec is a list of Ca2+ transient states (0/1 as a function of time)

PosCounts = zeros(NumPosBins);
ActCounts = zeros(NumActBins);

for i = 1:length(PosVec)
    PosCounts(PosVec(i)) = PosCounts(PosVec(i)) + 1;
    ActCounts(ActVec(i)) = ActCounts(ActVec(i)) + 1;
end

for i = 1:NumPosBins
    pPos(i) = PosCounts(i)/length(PosVec);
end

for i = 1:NumActBins
    pAct(i) = ActCounts(i)/length(ActVec);
end

I = 0;

for i = 1:NumPosBins
    for j = 1:NumActBins
        JP(i,j) = sum((PosVec == i) .* (ActVec == j))/length(ActVec);
        I = I + JP(i,j)*log2(JP(i,j)./pPos(i)./pAct(j));
    end
end
keyboard;
end

