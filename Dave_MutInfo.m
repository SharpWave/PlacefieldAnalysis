function [I] = Dave_MutInfo(PosVec,NumPosBins,ActVec,NumActBins)
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
        if ((pPos(i) > 0) & (pAct(j) > 0))
            JP(i,j) = sum((PosVec == i) .* (ActVec == j))/length(ActVec);
            if (JP(i,j) > 0)
            I = I + JP(i,j)*log2(JP(i,j)./pPos(i)./pAct(j));
            end
            if (isnan(I))
                keyboard;
            end
        end
    end
end

end

