function [I] = SkaggsCaMutInfo(TMap_unsmoothed,RunOccMap,NumRunFrames)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

I = 0;

NumSpatBins = length(TMap_unsmoothed(:));
RunOccP = RunOccMap./NumRunFrames;
TP = TMap_unsmoothed.*RunOccMap;
TP = TP./NumRunFrames;
TCounts = TP.*NumRunFrames;
FR = sum(TCounts(:))./NumRunFrames;

for i = 1:NumSpatBins
    if ((TP(i) > 0) && (RunOccP(i) > 0))
    I = I + RunOccP(i)*TMap_unsmoothed(i)*log2(TMap_unsmoothed(i)/FR);
    end
end

I = I./FR;

