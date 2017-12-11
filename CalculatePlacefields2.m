function [output_filename] = CalculatePlacefields2(varargin)
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
%       and 2nd half of session along with whole session maps
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
cmperbin = 0.25; % Dombeck uses 2.5 cm bins, Ziv uses 2x2 bins with 3.75 sigma gaussian smoothing
mouseradiuscm = 1.5;
use_mut_info = 0; % default
minspeed = 1; % cm/s, default
pos_align_file = '';
man_savename = [];
use_unaligned_data = 0; % default
alt_inputs = []; % default
HalfWindow = 0;
NumShuffles = 1000; % default for running bootstrapping in StrapIt below
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
        exclude_frames = varargin{j+1};
    end
    if strcmpi('exclude_frames_raw',varargin{j})
        exclude_frames_raw = varargin{j+1};
    end
    if strcmpi('rotate_to_std',varargin{j})
        rotate_to_std = varargin{j+1};
    end
    if strcmpi('name_append',varargin{j})
        name_append = ['_' varargin{j+1}];
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
    load(alt_inputs,'NeuronImage','NeuronPixelIdxList','PSAbool');
else
    load('FinalOutput.mat','NeuronImage','NeuronPixelIdxList','PSAbool');
    disp('Using FinalOutput.mat');
end
NeuronPixels = NeuronPixelIdxList;
FT = PSAbool;
NumNeurons = size(FT,1);
SR = 10,

try
    load('Pix2Cm.mat');
catch
    load('Pos.mat','ypos_interp','xpos_interp');
    figure(479);plot(ypos_interp,xpos_interp);axis image;
    FeatureCm = input('enter length in cm of known maze boundary ');
    FeaturePix = input('enter length in pixels of known maze boundary ');
    Pix2Cm = FeatureCm/FeaturePix;
    save Pix2Cm.mat Pix2Cm;
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
    
    [x,y,speed,FT,FToffset,FToffsetRear, aviFrame] = AlignImagingToTracking(Pix2Cm,FT,HalfWindow,SR);
    xmax = max(x); xmin = min(x);
    ymax = max(y); ymin = min(y);
    pos_align_use = 0;
end

Flength = length(x);

%% Adjust exclude_frames_raw if applicable

if ~isempty(exclude_frames_raw)
    % take raw/non-aligned frame inidices that are aligned with FT from
    % ProcOut.mat file and align to aligned position/FT data.
    exclude_aligned = exclude_frames_raw - FToffset + 2;
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
half_validonly = find(cumsum(ind_use) == round(sum(ind_use)/2),1,'last'); % Get halfway point of valid indices
ind_use_half{1}(1:half) = 1; %zeros(1,length(1:half)); % Get 1st half valid indices
ind_use_half{2}(half+1:length(ind_use)) = 1; %zeros(1,length(half+1:length(ind_use))); % get 2nd half valid indices
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

figure(13);plot(smspeed);minspeed,

runepochs = NP_FindSupraThresholdEpochs(smspeed,minspeed);
isrunning = smspeed >= minspeed;

t = (1:length(x))./SR;

figure(1);plot(t,smspeed);axis tight;xlabel('time (sec)');ylabel('speed cm/s');

% Set up binning and smoothers for place field analysis

% Dombeck used 2.5 cm bins - Ziv uses 2cm bins
xbuf = 0.15*(xmax-xmin);
xmax = round(xmax+xbuf);
xmin = round(xmin-xbuf);
ymax = round(ymax+xbuf);
ymin = round(ymin-xbuf);

Xrange = xmax-xmin;
Yrange = ymax-ymin;

NumXBins = ceil(Xrange/cmperbin);
NumYBins = ceil(Yrange/cmperbin);

Xedges = (0:NumXBins)*cmperbin+xmin;
Yedges = (0:NumYBins)*cmperbin+ymin;

figure(12);hold on;plot(x,y);title('animal trajectory');

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

Xbin(find(Xbin == 0)) = 1; % This is wrong and will send out of range values to 1 for some reason...not good!
Ybin(find(Ybin == 0)) = 1;

PositionVector = sub2ind([NumXBins,NumYBins],Xbin,Ybin);

NumPos = NumXBins * NumYBins;

NumNeurons = size(FT,1);
NumSamples = size(FT,2);
LastPV = 0;
Curr = 0;

% first, make NewFT and new PV
templengths = NP_FindSupraThresholdEpochs(isrunning,eps,0);
startruns = templengths(:,1);
endruns = templengths(:,2);

runidx = zeros(1,NumSamples);

runlengths = endruns-startruns+1;

mousedisk = fspecial('disk',round(mouseradiuscm/cmperbin));
mouselen = size(mousedisk,1);
mouseoff = round(median(1:size(mousedisk,1)));
NewPosMap = zeros(NumXBins,NumYBins);

for i = 1:length(startruns);
    for t = startruns(i):endruns(i)
        NewPosMap(Xbin(t)-mouseoff:Xbin(t)-mouseoff+mouselen-1,Ybin(t)-mouseoff:Ybin(t)-mouseoff+mouselen-1) =  NewPosMap(Xbin(t)-mouseoff:Xbin(t)-mouseoff+mouselen-1,Ybin(t)-mouseoff:Ybin(t)-mouseoff+mouselen-1)+mousedisk;
    end
end
figure(30);imagesc(NewPosMap);axis image;colorbar;
%

NumShuffles = 500;

totalrun = sum(runlengths);

for i = 1:length(runlengths)
    fakebase{i} = zeros(1,runlengths(i));
    binedge(i) = sum(runlengths(1:i))/totalrun;
end

binedge = [0,binedge];
largestfree = runlengths;

 for i = 1:NumNeurons 
    i/NumNeurons,
    
    
    f = FT(i,:);
    
    % calculate lengths of firing episodes
    trlengths = [];
    for j = 1:length(startruns)
        temp = NP_FindSupraThresholdEpochs(f(startruns(j):endruns(j)),eps,0);
        if(~isempty(temp))
            trlengths = [trlengths;(temp(:,2)-temp(:,1))+1];
        end
    end
    str = sort(trlengths,'descend');
    
    % start shuffling
    for j = 1:NumShuffles
        tempbase = fakebase;
        j/NumShuffles;
        largestfree = runlengths;
        numfree = runlengths;
        totalfree = totalrun;
        fakef = zeros(size(f));
        for k = 1:length(str)
            NotPlaced = true;
            while(NotPlaced)
                [~,BinToTry] = histc(rand,binedge);
                if (str(k) <= largestfree(BinToTry))
                    % we can put this transient into this run
                    
                    % find legal spots to stick this thing
                    tempgap = NP_FindSupraThresholdEpochs(~tempbase{BinToTry},eps,0);
                    gaplen = tempgap(:,2)-tempgap(:,1)+1;
                    PosPos = [];
                    for m = 1:length(gaplen)
                        if(gaplen(m) >= str(k))
                            PosPos = [PosPos,tempgap(m,1):tempgap(m,1)+gaplen(m)-str(k)];
                        end
                    end
                    rPos = ceil(rand*length(PosPos));
                    
                    % insert the new fake transient
                    tempbase{BinToTry}(PosPos(rPos):PosPos(rPos)+str(k)-1) = 1;
                    
                    % recalculate largestfree
                    tempgap = NP_FindSupraThresholdEpochs(~tempbase{BinToTry},eps,0);
                    if(~isempty(tempgap))
                        gaplen = tempgap(:,2)-tempgap(:,1)+1;
                    else
                        gaplen = 0;
                    end
                    
                    largestfree(BinToTry) = max(gaplen);
                    numfree(BinToTry) = numfree(BinToTry)-str(k);
                    totalfree = totalfree-str(k);
                    binedge = [];
                    for m = 1:length(largestfree)
                        binedge(m) = sum(numfree(1:m))/totalfree;
                    end
                    binedge = [0,binedge];
                    NotPlaced = false;
                end
            end
        end
        % assemble fakef
        for k = 1:length(tempbase)
            tempbase{k};
            fakef(startruns(k):endruns(k)) = fakef(startruns(k):endruns(k)) | tempbase{k};
        end
        fakePF = zeros(NumXBins,NumYBins);
        for k = 1:length(tempbase)
            for m = startruns(k):endruns(k)
                if fakef(m)
                    fakePF(Xbin(m)-mouseoff:Xbin(m)-mouseoff+mouselen-1,Ybin(m)-mouseoff:Ybin(m)-mouseoff+mouselen-1) =  fakePF(Xbin(m)-mouseoff:Xbin(m)-mouseoff+mouselen-1,Ybin(m)-mouseoff:Ybin(m)-mouseoff+mouselen-1)+mousedisk;
                end
            end
        end
        fakePF = fakePF./NewPosMap;
        PLcube(:,:,j) = fakePF;
        
    end
    % Calculate real place field
    realPL{i} = zeros(NumXBins,NumYBins);
    for k = 1:length(tempbase)
        for m = startruns(k):endruns(k)
            if FT(i,m)
                realPL{i}(Xbin(m)-mouseoff:Xbin(m)-mouseoff+mouselen-1,Ybin(m)-mouseoff:Ybin(m)-mouseoff+mouselen-1) =  realPL{i}(Xbin(m)-mouseoff:Xbin(m)-mouseoff+mouselen-1,Ybin(m)-mouseoff:Ybin(m)-mouseoff+mouselen-1)+mousedisk;
            end
        end
    end
    
    realPL{i} = realPL{i}./NewPosMap;
    PLpct{i} = zeros(NumXBins,NumYBins);
    
    for j = 1:NumXBins
        for k = 1:NumYBins
            svals = sort(PLcube(j,k,:));
            [~,idx] = min(abs(svals-realPL{i}(j,k)));
            PLpct{i}(j,k) = idx/NumShuffles;
            PLthr{i}(j,k) = realPL{i}(j,k) > mean(svals)+2*std(svals);
        end
    end
    %PLpct{i}(NewPosMap < 1) = 0;
    
    
    xb = 16; % this is a hack need to determine this auto
    sigPF{i} = realPL{i};
    sigPF{i}(PLthr{i} == 0) = 0;
    sigPF{i}(1:xb,:) = 0;
    sigPF{i}(end-xb+1:end,:) = 0;
    sigPF{i}(:,1:xb) = 0;
    sigPF{i}(:,end-xb+1:end) = 0;
    %temp(countPL{i} < 1.5) = 0;
%     imagesc(temp);colorbar;
%     hold on;
%     plot(Ybin,Xbin,'Color',[0.7 0.7 0.7])
%     a = find(FT(i,:).*isrunning);
%     plot(Ybin(a),Xbin(a),'ro');
%     a = find(FT(i,:).*~isrunning);
%     plot(Ybin(a),Xbin(a),'ko');
%     hold off;
%     axis image;
%     set(gcf,'Position',[20   116   634   854]);
%     figure(3);imagesc(PLpct{i});set(gca,'YDir','reverse');colorbar;axis image;set(gcf,'Position',[658   116   634   854]);
%     %figure(5);imagesc(countPL{i});set(gca,'YDir','reverse');colorbar;axis image;
%     pause;
end

save PlaceMaps2.mat sigPF NewPosMap realPL PLpct PLthr FT Xbin Ybin NumXBins NumYBins isrunning smspeed Xedges Yedges x y NeuronImage;




