function [] = PFstats(rot_to_std )
%PFstats(rot_to_std )
%   Calculate statistics on place-fields
%
% INPUTS
%
%   rot_to_sta: 0(default) uses data that has either not been aligned with
%   other data or not rotated such that local cues align, 1 - uses data
%   that has been rotated such that local cues align

if nargin == 0
    rot_to_std = 0;
end

% Load appropriate PlaceMaps file
if rot_to_std == 0
    load PlaceMaps.mat; % x y t xOutline yOutline speed minspeed FT TMap RunOccMap OccMap SpeedMap RunSpeedMap NeuronImage NeuronPixels cmperbin pval Xbin Ybin;
elseif rot_to_std == 1
   load PlaceMaps_rot_to_std.mat; 
end

% Which pixels are in place field?

% what should the threshold be?  peaks range from .1 to .35
% i.e. the BEST place cell is active in 35% of crossings?


NumNeurons = length(NeuronImage);
NumFrames = length(Xbin);

% some analysis using bwconncomp and regionprops

for i = 1:NumNeurons
    display(['calculating PF center for neuron ',int2str(i)]);
    peakval = max(TMap{i}(:));
    ThreshMap = TMap{i}.*(TMap{i} > peakval/2);
    BoolMap = ThreshMap > 0;
    b{i} = bwconncomp(BoolMap);
    r{i} = regionprops(b{i},'area','centroid');
end



for i = 1:NumNeurons
    display(['repackaging PF info for neuron ',int2str(i)])
    NumPF(i) = b{i}.NumObjects;
    PFpixels{i,1} = [];
    PFcentroid{i,1} = [];
    PFsize(i,1) = 0;
    MaxPF(i) = 0;
    for j = 1:NumPF(i)
        PFpixels{i,j} = b{i}.PixelIdxList{j};
        PFsize(i,j) = r{i}(j).Area;
        PFcentroid{i,j} = r{i}(j).Centroid;
    end
    [~,MaxPF(i)] = max(PFsize(i,:));
    % NK Question - should we calculate something similar for
    % firing/transient rate? e.g. if we have a cell with a very large field
    % that the cell only fires in once, but another, smaller field that the
    % cell fires reliably and on numerous passes, wouldn't it be good to
    % know that also?
end

% for every place field, figure out how many times the mouse passed through
% it

% convert Xbin and Ybin into a single number denoting where the mouse is
loc_index = sub2ind(size(TMap{1}),Xbin,Ybin);
for i = 1:NumNeurons
    display(['calculating PF visits for neuron ',int2str(i)])
    PFnumepochs(i,1) = 0;
    PFepochs{i,1} = [];
    for j = 1:NumPF(i)
        PixelBool = zeros(1,NumFrames);
        minp = min(PFpixels{i,j});
        maxp = max(PFpixels{i,j});
        for k = 1:NumFrames
            if (loc_index(k) < minp)
                PixelBool(k) = 0;
                continue;
            end
            if (loc_index(k) > maxp)
                PixelBool(k) = 0;
                continue;
            end
            PixelBool(k) = ismember(loc_index(k),PFpixels{i,j});
        end
        PFepochs{i,j} = NP_FindSupraThresholdEpochs(PixelBool,eps,0);
        PFnumepochs(i,j) = size(PFepochs{i,j},1);
    end
end

for i = 1:NumNeurons
    display(['calculating PF hits for neuron ',int2str(i)])

    for j = 1:NumPF(i)
        PFactive{i,j} = [];
        PFnumhits(i,j) = 0;
        PFpcthits(i,j) = 0;
        for k = 1:PFnumepochs(i,j)
            PFactive{i,j}(k) = sum(FT(i,PFepochs{i,j}(k,1):PFepochs{i,j}(k,2))) > 0;
        end
        if (~isempty(PFactive{i,j}))
            PFnumhits(i,j) = sum(PFactive{i,j});
            PFpcthits(i,j) = PFnumhits(i,j)/PFnumepochs(i,j);
        end
    end
end

if rot_to_std == 0
    save PFstats.mat PFpcthits PFnumhits PFactive PFnumepochs PFepochs MaxPF PFcentroid PFsize PFpixels -v7.3;
elseif rot_to_std == 1
    save PFstats_rot_to_std.mat PFpcthits PFnumhits PFactive PFnumepochs PFepochs MaxPF PFcentroid PFsize PFpixels -v7.3;
end

end





