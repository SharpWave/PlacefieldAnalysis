function [output_filename] = CalculatePlacefields(RoomStr,varargin)
% [] = [] = CalculatePlacefields(RoomStr,varargin)
% RoomStr, e.g. '201a'
%       Takes tracking data from Pos.mat (or Pos_align.mat, see batch_pos_align)
%       along with neural imaging calcium transient data from ProcOut.mat
%       in the working directory and calculates the transient heat map for each n
%       neuron in Proc.mat. and
%       determines if the field is statistically significant via
%       bootstrapping each placefield's entropy.
%
% varargins
%       -'progress_bar': 1 uses a progress bar in lieu of spam to screen
%       while running StrapIt (needs ProgressBar function written by Stefan
%       Doerr).  Default = 1;
%
%       -'exclude_frames': 1 x n array of frame numbers you wish to exclude from
%       PFA analysis. IMPORTANT - these frames correspond to position/FT
%       data that has already been aligned by running
%       AlignImagingToTracking.
%
%       -'exclude_frames_raw: same as exclude_frames but using indices from
%       raw FT data out of ProcOut.mat that has NOT been aligned to
%       position data (e.g. bad/dropped frames identified in Mosaic).  Can
%       be used in conjunction with 'exclude_frames' if you have both types
%       of frames you wish to exclude
%
%       -'rotate_to_std': 1 =  use position data that has been rotated back
%       such that all local cues are aligned (found in Pos_align_corr_std.mat). 
%       if using pre-aligned data (run batch_align_pos)
%       0 (default) = use data that has been aligned but NOT rotated back so that 
%       distal cues align (found in Pos_align.mat)
%
%       -'name_append': enter in the name you want to append to the
%       'PlaceMaps' file
%
%       -'cmperbin': centimeters per occupancy bin.  Default is 1.
%
%       -'calc_half': 0 = default. 1 = calculate TMap and pvalues for 1st
%       and 2nd half of session along with whole session maps. 2 = calc odd
%       versus even minutes.
%
%       -'mispeed': threshold for calculating placemaps.  Any values below
%           are not used. 1 cm/s = default.  
%
%       -'pos_align_file': use to load a Pos_align file that is not either
%       Pos_align.mat or Pos_align_std_corr.mat. must follow
%       'pos_align_file' with two arguments: 1) the name of the file to
%       load, and 2) a name to append to the Placefields file that will be
%       saved as output
%       
%       -'man_savename': use specified savename. 
%
%       -'use_unaligned_data': use Pos.mat file in lieu of Pos_align files
%       (default is to use aligned data!)
close all;

progress_bar = 1;
exclude_frames = [];
exclude_frames_raw = [];
rotate_to_std = 0;
name_append = '';
name_append2 = '';
calc_half = 0;
cmperbin = 1; % Dombeck uses 2.5 cm bins, Ziv uses 2x2 bins with 3.75 sigma gaussian smoothing
use_mut_info = 0; % default
minspeed = 1; % cm/s, default
pos_align_file = '';
man_savename = [];
use_unaligned_data = 0; % default
alt_inputs = []; % default
HalfWindow = 0; % default
NumShuffles = 500; % default for running bootstrapping in StrapIt below
for j = 1:length(varargin)
    if strcmpi('half_window',varargin{j})
        HalfWindow = varargin{j+1};
    end
    if strcmpi('progress_bar',varargin{j})
        progress_bar = varargin{j+1};
    end
    if strcmpi('man_savename',varargin{j})
        man_savename = varargin{j+1};
    end
    if strcmpi('alt_inputs',varargin{j})
        alt_inputs = varargin{j+1};
    end
    if strcmpi('exclude_frames',varargin{j})
        if ~isempty(exclude_frames)
            exclude_frames = varargin{j+1};
        end
    end
    if strcmpi('exclude_frames_raw',varargin{j})
        exclude_frames_raw = varargin{j+1};
    end
    if strcmpi('rotate_to_std',varargin{j})
        rotate_to_std = varargin{j+1};
    end
    if strcmpi('name_append',varargin{j})
        name_append = varargin{j+1};
    end
    if strcmpi('cmperbin',varargin{j})
       cmperbin = varargin{j+1};
    end
    if strcmpi('calc_half',varargin{j})
       calc_half = varargin{j+1};
    end
    if strcmpi(varargin{j},'use_mut_info')
        use_mut_info = varargin{j+1};
    end
    if strcmpi(varargin{j},'minspeed')
        minspeed = varargin{j+1};
    end
    if strcmpi(varargin{j},'pos_align_file')
        pos_align_file = varargin{j+1};
        name_append2 = varargin{j+2};
    end
    if strcmpi(varargin{j},'use_unaligned_data')
        use_unaligned_data = varargin{j+1};
    end
    if strcmpi(varargin{j},'NumShuffles')
        NumShuffles = varargin{j+1};
    end
end
name_append = [name_append name_append2];

if (~isempty(alt_inputs))
    display('using user-specified neuronal activity');
    load(alt_inputs,'NeuronImage','NeuronPixels','FT');
else
    load ('ProcOut.mat','NeuronImage','NeuronPixels','FT');
end

NumNeurons = size(FT,1);

SR = 20;
Pix2Cm = 0.15;

% Note that Pix2Cm should probably live in MakeMouseSession somewhere since
% there are actually a wide variety of these values.  Below is good enough
% for now, but the two-environment task in Aug/Sep 2015 will require a
% different value for sure
if (nargin == 0)
    Pix2Cm = 0.15;
    display('assuming room 201b');
    % factor for 201a is 0.0709
else
    if (strcmpi(RoomStr,'201a'))
        Pix2Cm = 0.0709;
        display('Room 201a');
    elseif strcmpi(RoomStr,'201b')
        Pix2Cm = 0.15;
        display('Room 201b');
    elseif strcmpi(RoomStr,'201a - 2015')
        Pix2Cm = 0.0874;
        display('Room 201a - 2015');
    end
end

for i = 1:NumNeurons
    temp = bwboundaries(NeuronImage{i});
    yOutline{i} = temp{1}(:,1);
    xOutline{i} = temp{1}(:,2);
end

try % Pull aligned data
    if ~isempty(pos_align_file) % Load alternate file with aligned position data if specified
        load(pos_align_file)
        disp(['Using position data that has been aligned to other like sessions in file "' pos_align_file '"'])
    elseif isempty(pos_align_file) && use_unaligned_data == 0  % load either Pos_align or Pos_align_std_corr
        if rotate_to_std == 0
            load Pos_align.mat
            disp('Using position data that has been aligned to other like sessions.')
        elseif rotate_to_std == 1
            load Pos_align_std_corr.mat
            disp('Using position data that has been aligned to other like sessions AND rotated so that local cues are aligned.')
        end
    elseif isempty(pos_align_file) && use_unaligned_data == 1
            disp('Loading Pos.mat for position data that has NOT been aligned to other sessions')
    end
    % Note that xmin, xmax, ymin, and ymax (used below) have been pulled from
    % Pos_align.mat.
    x = x_adj_cm;
    y = y_adj_cm;
    pos_align_use = 1;
    if ~exist('aviFrame','var')
        disp('Faking aviFrame variable - re-run AlignImagingToTracking to fix')
        aviFrame = 'not present - re-run AlignImagingToTracking to fix';
    end
    
    
catch % If no alignment has been performed, alert the user
    disp('Using position data that has NOT been aligned to other like sessions.')
    disp('NOT good for comparisons across sessions...run batch_align_pos for this.')

    [x,y,speed,FT,FToffset,FToffsetRear, aviFrame] = AlignImagingToTracking(Pix2Cm,FT,HalfWindow);
    xmax = max(x); xmin = min(x);
    ymax = max(y); ymin = min(y);
    pos_align_use = 0;
end

Flength = length(x);

%% Adjust exclude_frames_raw if applicable

if ~isempty(exclude_frames_raw)
    % take raw/non-aligned frame inidices that aligned with FT from
    % ProcOut.mat file and align to aligned position/FT data.
    load('Pos.mat','start_time');
    exclude_aligned = exclude_frames_raw - FToffset + round(start_time*SR); % Account for offsets between a) AVI start time, b) imaging start time, and c) Mouse on Maze time
    exclude_aligned = exclude_aligned(exclude_aligned > 0); % Get rid of any negative values (corresponding to times before the mouse was on the maze)
    exclude_frames = [exclude_frames, exclude_aligned]; % concatenate to exclude_frames 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get vector of valid frames to use

% If using aligned data, exclude any position data from outside the
% specified limits (e.g. if you only want to look at one part of an arena).
pos_ind_use = ones(1,Flength);
if pos_align_use == 1
    pos_ind_exclude = find(x < xmin | x > xmax | y < ymin | y > ymax);
    pos_ind_use(pos_ind_exclude) = zeros(1,length(pos_ind_exclude));
end
% Exclude frames specified in the session structure
ind_use = ones(1,Flength);
ind_use_half{1} = zeros(1,Flength);
ind_use_half{2} = zeros(1,Flength);
half = round(Flength/2);
% Check for edge case where exclude frames extend beyond the end of FT
% (usually due to adjustments in case of smoothing) - COULD BE A BUG HERE -
% DOES exclude_frames get messed up/ misaligned by half the smoothing
% window as a result of this???
if max(exclude_frames) > Flength
   temp = exclude_frames(exclude_frames <= Flength);
   exclude_frames = temp;
end
ind_use(exclude_frames) = zeros(1,length(exclude_frames)); % Send bad frames to zero
half_validonly = find(cumsum(ind_use) == round(sum(ind_use)/2),1,'last'); % Get halfway point of valid indices - this is NOT right yet, needs to get start and end of valid indices too

if calc_half == 1 % 1st half v 2nd half
    ind_use_half{1}(1:half) = 1; %zeros(1,length(1:half)); % Get 1st half valid indices
    ind_use_half{2}(half + 1:length(ind_use)) = 1; %zeros(1,length(half+1:length(ind_use))); % get 2nd half valid indices
elseif calc_half == 2 % odd v even minutes
    [odd_ind, even_ind] = odd_v_even_minutes(1:length(ind_use),SR);
    ind_use_half{1}(odd_ind) = 1;
    ind_use_half{2}(even_ind) = 1;
end
% Use only frames that are not excluded in the session structure AND are
% within the specified arena limits
frames_use_ind = ind_use & pos_ind_use;
frames_use_ind_half{1} = ind_use_half{1} & pos_ind_use; % NRK note: this should probably be changed to chop only the valid indices in half, not do an AND between the first half indices and valid indices
frames_use_ind_half{2} = ind_use_half{2} & pos_ind_use;
frames_use = find(frames_use_ind); 

%%%%%%%%%%%%%%%%%%%%%%%%

smspeed = convtrim(speed,ones(1,2*SR))./(2*SR);
% Adjust smspeed to include ONLY frames_use
temp = zeros(size(smspeed));
temp(frames_use_ind) = smspeed(frames_use_ind);
smspeed = temp;

runepochs = NP_FindSupraThresholdEpochs(smspeed,minspeed);
isrunning = smspeed >= minspeed;

t = (1:length(x))./SR;

figure(1);plot(t,smspeed);axis tight;xlabel('time (sec)');ylabel('speed cm/s');

% Set up binning and smoothers for place field analysis

% Dombeck used 2.5 cm bins - Ziv uses 2cm bins
Xrange = xmax-xmin;
Yrange = ymax-ymin;

NumXBins = ceil(Xrange/cmperbin);
NumYBins = ceil(Yrange/cmperbin);

Xedges = (0:NumXBins)*cmperbin+xmin;
Yedges = (0:NumYBins)*cmperbin+ymin;

figure(2);hold on;plot(x,y);title('animal trajectory');

% draw all of the edges
for i = 1:length(Xedges)
    z = line([Xedges(i) Xedges(i)],[Yedges(1) Yedges(end)]);
    set(z,'Color','r');
end

for i = 1:length(Yedges)
    z = line([Xedges(1) Xedges(end)],[Yedges(i) Yedges(i)]);
    set(z,'Color','r');
end

axis tight;

% Find out which bin the mouse was in at which timepoint

[~,Xbin] = histc(x,Xedges);
[~,Ybin] = histc(y,Yedges);

Xbin(find(Xbin == (NumXBins+1))) = NumXBins;
Ybin(find(Ybin == (NumYBins+1))) = NumYBins;

Xbin(find(Xbin == 0)) = 1;
Ybin(find(Ybin == 0)) = 1;

PositionVector = sub2ind([NumXBins,NumYBins],Xbin,Ybin);

RunOccMap = zeros(NumXBins,NumYBins); % # of samples in bin while running
RunOccMap_half{1} = zeros(NumXBins,NumYBins); %RunOccMap for 1st half of session
RunOccMap_half{2} = zeros(NumXBins,NumYBins); %RunOccMap for 2nd half of session
OccMap = zeros(NumXBins,NumYBins); % total # of samples in bin
OccMap_half{1} = zeros(NumXBins,NumYBins); % total # of samples in bin - 1st half
OccMap_half{2} = zeros(NumXBins,NumYBins); % total # of samples in bin - 2nd half
SpeedMap = zeros(NumXBins,NumYBins); % average speed in bin
RunSpeedMap = zeros(NumXBins,NumYBins); % average speed in bin while running

% Calculate Occupancy maps, both for all times and for times limited to 
% when the mouse was moving above minspeed
% OccMap and RunOccMap are in # of visits
first_half = frames_use(floor(length(frames_use))); % Get halfway point of session
for j = 1:length(frames_use)
    i = frames_use(j); % Grab next good frame to use
    if (isrunning(i))
        RunOccMap(Xbin(i),Ybin(i)) = RunOccMap(Xbin(i),Ybin(i))+1;
        if (i ~= Flength)
            RunSpeedMap(Xbin(i),Ybin(i)) = RunSpeedMap(Xbin(i),Ybin(i))+smspeed(i);
        end
    end
    
    OccMap(Xbin(i),Ybin(i)) = OccMap(Xbin(i),Ybin(i))+1;
    if (i ~= Flength)
        SpeedMap(Xbin(i),Ybin(i)) = SpeedMap(Xbin(i),Ybin(i))+smspeed(i);
    end
    
    
    if i == first_half
        RunOccMap_half{1} = RunOccMap; % Save 1st half RunOccMap
        OccMap_half{i} = OccMap;
    elseif i == length(frames_use)
        RunOccMap_half{2} = RunOccMap - RunOccMap_half{1}; % Save 2nd half RunOccMap
        OccMap_half{2} = OccMap - OccMap_half{1};
    end

end
SpeedMap = SpeedMap./OccMap;
RunSpeedMap = RunSpeedMap./RunOccMap;

p = ProgressBar(NumNeurons);

for i = 1:NumNeurons
  [TMap{i}, TMap_gauss{i}, TMap_unsmoothed{i}] = calcmapdec(FT(i,:), ...
      RunOccMap, Xbin, Ybin, isrunning & frames_use_ind, cmperbin);
  [pval(i), pvalI(i),SpatialH(i)] = StrapIt(FT(i,:), RunOccMap, Xbin, Ybin, cmperbin, runepochs, isrunning & frames_use_ind,...
      0, 'suppress_output', progress_bar,'use_mut_info',use_mut_info,'NumShuffles',NumShuffles);
  if calc_half == 1 || calc_half == 2 % Calculate half-session TMaps and p-values
      for j = 1:2
          [TMap_half(j).Tmap{i}, TMap_half(j).TMap_gauss{i}, TMap_half(j).TMap_unsmoothed{i}] = ...
              calcmapdec(FT(i,:), RunOccMap, Xbin, Ybin, isrunning & frames_use_ind_half{j}, cmperbin);
          [pval_half{j}.pval(i), pval_half{j}.pvalI(i),pval_half{j}.SpatialH(i)] = StrapIt(FT(i,:), RunOccMap, Xbin, Ybin, cmperbin, runepochs, isrunning & frames_use_ind_half{j},...
              0, 'suppress_output', progress_bar,'use_mut_info',use_mut_info,'NumShuffles',NumShuffles);
      end
  else
      TMap_half = [];
      pval_half = [];
  end
  SpatialI(i) = Dave_MutInfo(PositionVector(find(isrunning)),NumXBins*NumYBins,FT(i,find(isrunning))+1,2);
  
  if progress_bar == 1
     p.progress; 
  end
end
p.stop;


%PFreview(FT,TMap,t,x,y,pval,ip,find(pval > 0.95)) this finds all of the
%decent placefields
if isempty(man_savename)
    if rotate_to_std == 0
        save_name = ['PlaceMaps' name_append '.mat'] ;
    elseif rotate_to_std == 1
        save_name = ['PlaceMaps_rot_to_std' name_append '.mat'];
    end
else
    save_name = man_savename;
end

%%% NRK - save 1st and 2nd half stuff here
% save PlaceMaps.mat x y t xOutline yOutline speed minspeed FT TMap RunOccMap OccMap SpeedMap RunSpeedMap NeuronImage NeuronPixels cmperbin pval Xbin Ybin FToffset FToffsetRear isrunning cmperbin Xedges Yedges; 
save(save_name,'x', 'y', 't', 'xOutline', 'yOutline', 'speed','minspeed', ...
    'FT', 'TMap','TMap_gauss', 'TMap_unsmoothed', 'RunOccMap', 'OccMap', ...
    'SpeedMap', 'RunSpeedMap', 'NeuronImage', 'NeuronPixels',...
    'cmperbin', 'pval', 'Xbin', 'Ybin', 'FToffset', 'FToffsetRear', 'isrunning',...
    'Xedges', 'Yedges','exclude_frames','aviFrame','TMap_half','pval_half',...
    'pvalI','SpatialH','Pix2Cm','SpatialI','-v7.3'); 

output_filename = save_name;

return;


