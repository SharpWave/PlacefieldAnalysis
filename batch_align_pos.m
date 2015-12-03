function [ ] = batch_align_pos(base_struct, reg_struct, varargin)
% batch_align_pos(base_struct, reg_struct, varargin)
%
% Aligns position data so that every session has the same bounds on the
% occupancy map and we can easily do correlations and other comparisons
% between firing across sessions.  Everything gets scaled to the
% trajectory/occupancy from the base session. IMPORTANT NOTE: running this
% WILL overwrite any previous Pos_align.mat files, so be careful (though
% not too careful since using this function is fairly quick)
%
% INPUTS: 
%   base_struct & reg_struct:
%        mirror MD from MakeMouseSessionList, but must include at least
%       .Animal, .Date, .Session, AND .Room fields (and .Notes if you wish
%       to perform auto-rotation of the arena back to the standard
%       configuration)
%
% OPTIONAL INPUTS (specify as batch_align_pos(...,'manual_rot_overwrite,1,...):
%
%   manual_rot_overwrite: default = 1, will prompt you to
%       perform the rotation CORRECTION (e.g., re-aligning a slightly skewed
%       session) for each session, 0 will use pre-existing rotation
%       data in rotated.mat file in the working directory if it exists
%
%   ratio_use: ratio of the data to use for alignment - if 0.95, then
%       data is scaled so that the middle 95% of it in each session aligns with
%       the middle 95% in other sessions.  default = 0.95.
%
%   auto_rotate_to_std: set to 1 if you want to automatically rotate the
%       data back to the standard orientation if you have rotated the arena
%       as part of a control for distal v local cues.  Set to 0 if you want
%       to leave it alone (default).  IMPORTANT NOTES: 1) for this to work, the rotation of
%       the arena must be in the .Notes field of base_struct and
%       reg_struct. 2) 
%
%   manual_limits: set to 1 to manually draw the limits of the arena you
%       want to scale (default = 0).  Use this if, for example, you want to 
%       only consider the left part of a large arena that is joined to another.  
%       Input is either a single logical value, in which
%       case it applies to all sessions, or you can input a logical array where
%       the 1st entry corresponds to the base session, and the subsequent
%       entries correspond to all the registered sessions (e.g. [0 1 0 1] =
%       manually draw limits for the 1st and 3rd registered sessions but
%       not the base session or 2nd registered session)
%
% OUTPUTS (saved in Pos_align.mat in working directory, or Pos_align_std_corr.mat
%          if you choose to auto-rotate back):
%
%   'x_adj_cm','y_adj_cm': x and y positions converted to centimeters and
%   scaled/aligned to the base session such that all sessions align
%
%   'xmin','xmax','ymin','ymax': min and max position data for ALL sessions
%   that can be sent to CalculatePlacefields such that the occupancy maps,
%   heat maps, etc. are all identical in size
%
%   'speed','FT','FToffset','FToffsetRear': calculated from  
%    
%   'base_struct':
%   'sessions_included':
%
% NOTE: this may not work well for the 2 environment experiment since it
% does not account for any fish-eye distortions of the maze...should be
% good for most comparisons between the same mazes, however

close all

%% Parameters/default values
manual_rot_overwrite = 1;
ratio_use = 0.95; % ratio of the data to use for alignment - if 0.95, then
% data is scaled so that the middle 95% of it in each session aligns with
% the middle 95% in other sessions
auto_rotate_to_std = 0;
manual_limits = zeros(1,length(reg_struct) + 1);
xmin = 10; % where you want to set the minimum x value
ymin = 20; % where you want to set the minimum y-value

%% 0: Get varargins

for j = 1:length(varargin)
   if strcmpi(varargin{j},'manual_rot_overwrite')
       manual_rot_overwrite = varargin{j+1};
   end
   if strcmpi(varargin{j},'ratio_use')
       ratio_use = varargin{j+1};
   end
   if strcmpi(varargin{j},'auto_rotate_to_std')
       auto_rotate_to_std = varargin{j+1};
   end
   if strcmpi(varargin{j},'manual_limits')
      manual_limits = varargin{j+1}; 
   end
end

%% 1: Load all sessions, and align to imaging data

% Dump everything into one structure for future ease
sesh(1) = base_struct;
sesh(2:length(reg_struct) + 1) = reg_struct;

currdir = cd;
for j = 1: length(sesh)
    ChangeDirectory(sesh(j).Animal, sesh(j).Date ,sesh(j).Session);
    if ~isempty(regexpi(sesh(1).Room,'201b'))
        Pix2Cm = 0.15;
        disp(['Using 0.15 for Pix2Cm for ' sesh(j).Date ' Session ' num2str(sesh(j).Session)])
    else
        Pix2Cm = [];
        disp('Need room to get Pix2Cm')
    end
    load('ProcOut.mat', 'FT')
    % Align tracking and imaging
    [x,y,speed,FT,FToffset,FToffsetRear, aviFrame] = AlignImagingToTracking(Pix2Cm,FT);
    
    % Auto-rotate back to standard configuration if indicated
    if auto_rotate_to_std == 1
        rot_corr = get_rot_from_db(sesh(j));
        [x, y] = rotate_arena(x,y,rot_corr);
    end
    sesh(j).x = x;
    sesh(j).y = y;
    sesh(j).FT = FT;
    sesh(j).speed = speed;
    sesh(j).FToffset = FToffset;
    sesh(j).FToffsetRear = FToffsetRear;
    % Fix day-to-day mis-alignments in rotation of the maze
    [~,rot_x,rot_y, rot_ang] = sections(x,y,0,'manual_rot_overwrite',manual_rot_overwrite);
    sesh(j).rot_x = rot_x;
    sesh(j).rot_y = rot_y;
    sesh(j).rot_ang = rot_ang;
    sesh(j).aviFrame = aviFrame;
    
end

% keyboard

%% 2: Align position data for each session to the base session by using the 95% occupancy limits, save as Pos_align.mat
% Include base session in Pos_align for future reference

% Add in opportunity to manually select data limits to use here... if the
% flag you set for it is 1 - may want to do this by looking at the Notes
% section in MD

for j = 1:length(sesh)
    
    if manual_limits(j) == 0
        x_for_limits = sesh(j).rot_x;
        y_for_limits = sesh(j).rot_y;
        sesh(j).ind_keep = logical(ones(1,length(sesh(j).rot_x)));
    elseif manual_limits(j) == 1
        [x_for_limits, y_for_limits, sesh(j).ind_keep] = draw_manual_limits(...
            sesh(j).rot_x, sesh(j).rot_y);
    end
    % Get ecdfs of all x and y points
    [sesh(j).e_fx, sesh(j).e_x] = ecdf(x_for_limits);
    [sesh(j).e_fy, sesh(j).e_y] = ecdf(y_for_limits);
    % Find limits that correspond to ratio_use (e.g. if ratio_use = 0.95,
    % look for the x value that corresponds to 0.025 and 0.975)
    xbound{j}(1) = sesh(j).e_x(findclosest((1-ratio_use)/2,sesh(j).e_fx));
    xbound{j}(2) = sesh(j).e_x(findclosest(1 - (1-ratio_use)/2,sesh(j).e_fx));
    ybound{j}(1) = sesh(j).e_y(findclosest((1-ratio_use)/2,sesh(j).e_fy));
    ybound{j}(2) = sesh(j).e_y(findclosest(1 - (1-ratio_use)/2,sesh(j).e_fy));
    % Calculate the span and get the ratio to the base span
    span_x(j) = xbound{j}(2) - xbound{j}(1);
    span_y(j) = ybound{j}(2) - ybound{j}(1);
    if j == 1
        span_x_ratio = 1;
        span_y_ratio = 1;
    elseif j > 1
        span_x_ratio = span_x(j)/span_x(1);
        span_y_ratio = span_y(j)/span_y(1);
    end
    
    % Linearly adjust all the coordinates to match - use all position data!
    sesh(j).x_adj = (sesh(j).rot_x - xbound{j}(1))/span_x_ratio + xmin;
    sesh(j).y_adj = (sesh(j).rot_y - ybound{j}(1))/span_y_ratio + ymin;
    
end
%% 4: Concatenate ALL position data into one X and one Y vector, and get Xedges and Yedges based on this

x_all = [];
y_all = [];
for j = 1:length(sesh)
    x_all = [x_all sesh(j).x_adj(sesh(j).ind_keep)]; % Only use data within the limits you drew for this!
    y_all = [y_all sesh(j).y_adj(sesh(j).ind_keep)]; % Only use data within the limits you drew for this!
end

%% 5: Get xmin, xmax, ymin, and ymax

xmax = max(x_all);
xmin = min(x_all);
ymax = max(y_all);
ymin = min(y_all);

%% 6: Save Xedges, Yedges in base session for future reference along with all sessions aligned to it.
% Also save adjusted position data for future use...

sessions_included(1) = base_struct;
sessions_included(2:length(reg_struct) + 1) = reg_struct;

for j = 1:length(sesh)
    x_adj_cm = sesh(j).x_adj;
    y_adj_cm = sesh(j).y_adj;
    speed = sesh(j).speed;
    FT = sesh(j).FT;
    FToffset = sesh(j).FToffset;
    FToffsetRear = sesh(j).FToffsetRear;
    aviFrame = sesh(j).aviFrame;
    if auto_rotate_to_std == 0
    save(fullfile(sesh(j).Location,'\Pos_align.mat'),'x_adj_cm','y_adj_cm',...
        'xmin','xmax','ymin','ymax', 'speed', 'FT', 'FToffset', ...
        'FToffsetRear', 'aviFrame', 'base_struct','sessions_included','auto_rotate_to_std');
    elseif auto_rotate_to_std == 1
        % finish here - save as a different filename?
        save(fullfile(sesh(j).Location,'\Pos_align_std_corr.mat'),'x_adj_cm','y_adj_cm',...
        'xmin','xmax','ymin','ymax', 'speed', 'FT', 'FToffset', ...
        'FToffsetRear','aviFrame', 'base_struct', 'sessions_included', 'auto_rotate_to_std');
    end
end

%% 7: Plot everything as a check
figure(100);
for j = 1:length(sesh)
    % Plot on an individual subplot
    subplot_auto(length(sesh) + 1,j+1);
    plot(sesh(j).x_adj,sesh(j).y_adj);
    xlim([xmin xmax]); ylim([ymin ymax])
    title(['Session ' num2str(j)])
    % Plot everything on top of the other
    subplot_auto(length(sesh) + 1, 1);
    hold on
    plot(sesh(j).x_adj, sesh(j).y_adj);
    xlim([xmin xmax]); ylim([ymin ymax])
    hold off
    title('All Sessions')
end

end

