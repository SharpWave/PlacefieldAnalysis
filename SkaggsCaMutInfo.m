function [I] = SkaggsCaMutInfo(TMap_unsmoothed,RunOccMap,NumRunFrames,cmperbin)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

I = 0;
gauss_std = 2.5; % cm
gauss_std = gauss_std/cmperbin;
sm_gauss = fspecial('gaussian',[round(8*gauss_std,0), round(8*gauss_std,0)],gauss_std);
RunOccMap = imfilter(RunOccMap,sm_gauss);
RunOccP = RunOccMap./sum(RunOccMap(:));



NumSpatBins = length(TMap_unsmoothed(:));
RunOccP = RunOccMap./NumRunFrames;
TP = TMap_unsmoothed.*RunOccMap;
TP = TP./NumRunFrames;
TCounts = TP.*NumRunFrames;
FR = sum(TCounts(:))./NumRunFrames;

for i = 1:NumSpatBins
    if (TMap_unsmoothed(i) > 0)
        if (RunOccP(i) > 0)
    I = I + RunOccP(i)*TMap_unsmoothed(i)*log2(TMap_unsmoothed(i)/FR);
        end
    end
end

I = I/FR;



