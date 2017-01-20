function [TMap, TMap_gauss, TMap_unsmoothed] = calcmapdec(trace,occmap,Xbin,Ybin,goodsamples, cmperbin, varargin)
% trace is boolean vector of cell on/off
% occmap is the occupancy map
% Xbin is the x coordinate of the mouse
% Ybin is the y coordinate of the mouse
% goodsamples is a boolean, 0 if the sample doesn't get used
% IMPORTANT the 'occmap' parameter and the 'goodsamples' parameter need to be
% consistent
%
% varargins:
%   'gauss_std': standard deviation in centimenters to use for gaussian 
%    smoothing filter. Default = 2.5cm.
%
%   'disk_rad': radius of disk smoothing filter to use in centimeters.
%   Default = 4 cm.
%
%   'plot_filters': default = 0.  1 plots out both the disk and gaussian
%   filtered TMaps and calcium "spike" trains for troubleshooting.

% keyboard

%% Smoothing variables and others
disk_rad = 4; % cm
gauss_std = 2.5; % cm
plot_filters = 0; % default - don't plot anything
calc_disk_only = false; % default - calculcates both disk and gaussian filtered map
for j = 1:length(varargin)
   if strcmpi('disk_rad',varargin{j})
       disk_rad = varargin{j+1};
   end
   if strcmpi('gauss_std',varargin{j})
       gauss_std = varargin{j+1};
   end
   if strcmpi('plot_filters',varargin{j})
       plot_filters = varargin{j+1};
   end
   if strcmpi('calc_disk_only',varargin{j})
       calc_disk_only = varargin{j+1};
   end
end

% Convert from centimeters into bins
disk_rad = disk_rad/cmperbin;
gauss_std = gauss_std/cmperbin;

%% Calculate Occupancy and "firing rate" and smooth

TMap = zeros(size(occmap));
Flength = length(trace);
for j = 1:Flength
    if(goodsamples(j))
        TMap(Xbin(j),Ybin(j)) = TMap(Xbin(j),Ybin(j)) + (trace(j) > 0);
    end
end

TMap = TMap./occmap; % Get Transient rate map
TMap(isnan(TMap)) = 0; % Set NaNs to zero
TMap_unsmoothed = TMap; % Save unsmoothed TMap
TMap2 = TMap; % Save copy for gaussian smoothing
Tsum = sum(TMap(:));

% Filter specifications
if ~calc_disk_only
    sm_gauss = fspecial('gaussian',[round(8*gauss_std,0), round(8*gauss_std,0)],gauss_std);
end
sm = fspecial('disk',disk_rad);

%%%  FOR REFERENCE %%%
% Used with 0.25 cm/bin
% sm_gauss = fspecial('gaussian',[100 100],7.5);
% sm = fspecial('disk',15);
%
% Used with 2 cm/bin
% sm_gauss = fspecial('gaussian',[10 10],1.75);
% sm = fspecial('disk',2);
%
% Used with 1 cm/bin
% sm_gauss = fspecial('gaussian',[20 20],2.5-3.5);
% sm = fspecial('disk',4);

%%% END REFERENCE %%%

TMap = imfilter(TMap,sm); % Apply disk filter
TMap = TMap.*Tsum./sum(TMap(:)); % keep sum the same
if calc_disk_only
    TMap_gauss = nan(size(TMap));
elseif ~calc_disk_only
    TMap_gauss = imfilter(TMap2,sm_gauss); % Apply gaussian filter
    TMap_gauss = TMap_gauss.*Tsum./sum(TMap_gauss(:)); % keep sum the same
end

%% Plotting for troubleshooting

if plot_filters == 1
    figure; 
    % Plot each type of filter
    subplot(2,2,1); imagesc(TMap); subplot(2,2,2); imagesc(TMap_gauss); colormap jet
    subplot(2,2,3); % Plot calcium transient spike train, with spikes occuring above the velocity threshold in red stars
    plot(trace); hold on; 
    plot(find(goodsamples & trace > 0), trace(goodsamples & trace > 0),'r*')
end

end