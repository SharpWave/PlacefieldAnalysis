function [x,y,speed,FT,FToffset,FToffsetRear,aviFrame,time_interp,nframesinserted] = AlignImagingToTracking(Pix2Cm,FT,HalfWindow, varargin)
% [x,y,speed,FT,FToffset,FToffsetRear,aviFrame,time_interp,nframesinserted] = AlignImagingToTracking(Pix2Cm,FT,HalfWindow, varargin)
%   Aligns imaging and tracking data.
%
%   IMPORTANT NOTE: You MUST be in the directory where Pos.mat lives for
%   this session for this to work.  Need to fix this somehow in the future,
%   but for now it is what it is...currently attempting to do with
%   basedir where you can add in the directory where Pos.mat is

p = inputParser;
p.addRequired('Pix2CM', @isnumeric);
p.addRequired('FT', @(a) isnumeric(a) || islogical(a));
p.addOptional('HalfWindow',0, @(a) round(a) == a && a >= 0);
p.addParameter('suppress_output',false, @islogical);
p.addParameter('basedir', pwd, @ischar);
p.parse(Pix2Cm, FT, HalfWindow, varargin{:});

HalfWindow = p.Results.HalfWindow;
suppress_output = p.Results.suppress_output;
basedir = p.Results.basedir;

SR = 20;

try 
    load(fullfile(basedir,'Pos.mat'))
    x = xpos_interp;
    y = ypos_interp;
    if ~exist('nframesinserted','var') % Backward's compatibility fix/notification
        if ~suppress_output
            disp('IMPORTANT - Cineplex dropped frame analysis not found in Pos.mat.')
            disp('Check with PreProcessMousePosition_auto to ensure your data is properly aligned!')
        end
        nframesinserted = nan;
    end
catch
    vidFile = dir('*.DVT');
    [xpos_interp,ypos_interp,start_time,MoMtime,time_interp,AVItime_interp,nframesinserted]...
        = PreProcessMousePosition_auto(vidFile.name);
end
x = xpos_interp;
y = ypos_interp;


x = x.*Pix2Cm;
y = y.*Pix2Cm;
dx = diff(x);
dy = diff(y);
speed = hypot(dx,dy)*SR;

% align starts of Fluorescence video and Plexon tracking to the arrival of
% the mouse of the maze
fTime = (1:size(FT,2))/SR;
fStart = findclosest(MoMtime,fTime);
FT = FT(:,fStart:end);

plexTime = (0:length(x)-1)/SR+start_time;
pStart = findclosest(MoMtime,plexTime);
x = x(pStart:end);
y = y(pStart:end);
time_interp = time_interp(pStart:end);          %time_interp starts when mouse is on maze. 

if exist('AVItime_interp','var')
    aviFrame = AVItime_interp(pStart:end);      %aviFrame starts from beginning of AVI. 
else
    aviFrame = 1:size(FT,2);
    if ~suppress_output
        disp('aviFrame not found in Pos.mat.  Faking for now')
    end
end

speed = speed(pStart:end);
plexTime = (1:length(x))/SR;

Flength = size(FT,2);
FToffsetRear = 0;

% if Inscopix or Plexon is longer than the other, chop
if (length(plexTime) <= Flength)
    % Chop the FL
    FT  = FT(:,1:length(plexTime));
    FToffsetRear = length(plexTime) - Flength;
    Flength = length(plexTime);
    
else
    speed = speed(1:Flength);
    x = x(1:Flength);
    y = y(1:Flength);
    aviFrame = aviFrame(1:Flength);
    time_interp = time_interp(1:Flength);
end

if (length(speed) < length(x))
    speed(end+1) = 0;
end

speed(1:100) = 0; % a hack, but otherwise screwy things happen

%%%%%%%%% adjust by half the movie smoothing window
%HalfWindow = 10;

% shift position and speed right
x = [zeros(1,HalfWindow),x(1:end-HalfWindow)];
y = [zeros(1,HalfWindow),y(1:end-HalfWindow)];
aviFrame = [zeros(1,HalfWindow),aviFrame(1:end-HalfWindow)];
speed = [zeros(1,HalfWindow),speed(1:end-HalfWindow)];
time_interp = [zeros(1,HalfWindow),time_interp(1:end-HalfWindow)];

% chop the first HalfWindow
x = x(HalfWindow+1:end);
y = y(HalfWindow+1:end);
aviFrame = aviFrame(HalfWindow+1:end);
speed = speed(HalfWindow+1:end);
FT = FT(:,HalfWindow+1:end);
time_interp = time_interp(HalfWindow+1:end);

FToffset = fStart + HalfWindow + 1;

end