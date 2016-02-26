function [] = PFstats(rot_to_std, varargin)
%PFstats(rot_to_std, varargin )
%   Calculate statistics on place-fields
%
% INPUTS
%
%   rot_to_std: 0(default) uses data that has either not been aligned with
%   other data or not rotated such that local cues align, 1 - uses data
%   that has been rotated such that local cues align
%
% varargins
%   -'alt_file_use': use to load a Pos_align file that is not either
%       Pos_align.mat or Pos_align_std_corr.mat. must follow
%       'alt_file_use' with  the name of the file to
%       load and the text you want to append onto the end of the PFstats
%       output file
%
%   -'tmap_thresh_denom':  used to threshold the TMap to define your place-field
%       extents by taking the max TMap value, dividing by tmap_thresh_denom, and
%       defining anything above this value as in-field (default = 2);

if nargin == 0
    rot_to_std = 0;
end

alt_file = 0;
progress_bar = 1; % default
name_append = '';
tmap_thresh_denom = 2; % default
for j = 1:length(varargin)
    if strcmpi(varargin{j},'alt_file_use')
        alt_file_use = varargin{j+1};
        name_append = varargin{j+2};
        if ~isempty(alt_file_use)
            alt_file = 1;
        end
    end
   if strcmpi('progress_bar',varargin{j})
            progress_bar = varargin{j+1};
   end
   if strcmpi('tmap_thresh_denom',varargin{j})
       tmap_thresh_denom = varargin{j+1};
   end
end

% Load appropriate PlaceMaps file
if alt_file == 1
    load(alt_file_use)
elseif alt_file == 0
    if rot_to_std == 0
        load PlaceMaps.mat; % x y t xOutline yOutline speed minspeed FT TMap RunOccMap OccMap SpeedMap RunSpeedMap NeuronImage NeuronPixels cmperbin pval Xbin Ybin;
    elseif rot_to_std == 1
        load PlaceMaps_rot_to_std.mat;
    end
end

% Which pixels are in place field?

% what should the threshold be?  peaks range from .1 to .35
% i.e. the BEST place cell is active in 35% of crossings?


NumNeurons = length(NeuronImage)
NumFrames = length(Xbin);

% some analysis using bwconncomp and regionprops
disp('Calculating PF centers for all neurons')
if progress_bar == 1
    p = ProgressBar(NumNeurons);
end
for i = 1:NumNeurons
    if progress_bar == 1
        p.progress;
    else
        display(['calculating PF center for neuron ',int2str(i)]);
    end
    peakval = max(TMap{i}(:));
    ThreshMap = TMap{i}.*(TMap{i} > peakval/tmap_thresh_denom);
    BoolMap = ThreshMap > 0;
    b{i} = bwconncomp(BoolMap);
    r{i} = regionprops(b{i},'area','centroid');
    
end
if progress_bar == 1
    p.stop;
end

disp('repackaging PF info for all neurons')
if progress_bar == 1
    p = ProgressBar(NumNeurons);
end
for i = 1:NumNeurons
    if progress_bar == 1
        p.progress;
    else
        display(['repackaging PF info for neuron ',int2str(i)])
    end
    
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
if progress_bar == 1
    p.stop;
end

% for every place field, figure out how many times the mouse passed through
% it

% convert Xbin and Ybin into a single number denoting where the mouse is
loc_index = sub2ind(size(TMap{1}),Xbin,Ybin);

disp('Calculating PF visits for all neurons')
if progress_bar == 1
    p = ProgressBar(NumNeurons);
end
for i = 1:NumNeurons
    if progress_bar == 1
        p.progress;
    else
        display(['calculating PF visits for neuron ',int2str(i)])
    end
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
if progress_bar == 1
    p.stop;
end

disp('Calculating PF hits for all neurons')
if progress_bar == 1
    p = ProgressBar(NumNeurons);
end
for i = 1:NumNeurons
    
    if progress_bar == 1
        p.progress;
    else
        display(['calculating PF hits for neuron ',int2str(i)])
    end

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
if progress_bar == 1
    p.stop;
end

if rot_to_std == 0
    save_name = ['PFstats' name_append '.mat'] ;
elseif rot_to_std == 1
    save_name = ['PFstats_rot_to_std' name_append '.mat'];
end

save(save_name, 'PFpcthits', 'PFnumhits', 'PFactive', 'PFnumepochs', 'PFepochs',...
    'MaxPF', 'PFcentroid', 'PFsize', 'PFpixels', 'tmap_thresh_denom', '-v7.3');

end





