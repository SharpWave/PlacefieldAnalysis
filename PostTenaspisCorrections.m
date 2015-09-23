function [] = PostTenaspisCorrections()
% uses extracted traces to correct the activation matrix produced by
% Tenaspis.

load ProcOut.mat;
load DumbTraces.mat;

NumNeurons = size(FT,1);

PFT = zeros(size(FT));

for i = 1:NumNeurons
    TRtimes = NP_FindSupraThresholdEpochs(FT(i,:),eps);
    InThreshes = Dtrace(i,(TRtimes(:,1)));
    AvgThresh(i) = mean(InThreshes);
    PFT(i,:) = Dtrace(i,:) > AvgThresh(i);
end

save Post_FT.mat PFT AvgThresh;

