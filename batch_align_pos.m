function [ ] = batch_align_pos(base_struct, reg_struct)
% Aligns position data so that every session has the same bounds on the
% occupancy map and we can easily do correlations and other comparisons
% between firing across sessions.  Everything gets scaled to the
% trajectory/occupancy from the base session.
%
% INPUTS: mirror MD from MakeMouseSessionList, but must include at least
% .Animal, .Date, .Session, AND .Room fields
%
% NOTE: this will not work well for the 2 environment experiment since it
% does not account for any fish-eye distortions of the maze...should be
% good for most comparisons between the same mazes, however

%% Parameters
ratio_use = 0.95; % ratio of the data to use for alignment - if 0.95, then
% data is scaled so that the middle 95% of it in each session aligns with
% the middle 95% in other sessions
xmin = 10; % where you want to set the minimum x value
ymin = 10; % where you want to set the minimum y-value

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
    [x,y,speed,FT,FToffset,FToffsetRear] = AlignImagingToTracking(Pix2Cm,FT);
    sesh(j).x = x;
    sesh(j).y = y;
    sesh(j).FT = FT;
    sesh(j).speed = speed;
    sesh(j).FToffset = FToffset;
    sesh(j).FToffsetRear = FToffsetRear;
    % Fix day-to-day mis-alignments in rotation of the maze
    [~,rot_x,rot_y, rot_ang] = sections(x,y,0);
    sesh(j).rot_x = rot_x;
    sesh(j).rot_y = rot_y;
    sesh(j).rot_ang = rot_ang;
    
end

keyboard

%% 2: Align position data for each session to the base session by using the 95% occupancy limits, save as Pos_align.mat
% Include base session in Pos_align for future reference

for j = 1:length(sesh)
    % Get ecdfs of all x and y points
    [sesh(j).e_fx, sesh(j).e_x] = ecdf(sesh(j).rot_x);
    [sesh(j).e_fy, sesh(j).e_y] = ecdf(sesh(j).rot_y);
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
    
    % Linearly adjust all the coordinates to match
    sesh(j).x_adj = (sesh(j).rot_x - xbound{j}(1))*span_x_ratio + xmin;
    sesh(j).y_adj = (sesh(j).rot_y - ybound{j}(1))*span_y_ratio + ymin;
    
end
%% 4: Concatenate ALL position data into one X and one Y vector, and get Xedges and Yedges based on this

x_all = [];
y_all = [];
for j = 1:length(sesh)
    x_all = [x_all sesh(j).x_adj];
    y_all = [y_all sesh(j).y_adj];
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
    FT = sesh(j).FT;
    FTalign = sesh(j).FTalign;
    FToffset = sesh(j).FToffset;
    FToffsetRear = sesh(j).FToffsetRear;
    save(fullfile(sesh(j).Location,'\Pos_align.mat'),'x_adj_cm','y_adj_cm',...
        'xmin','xmax','ymin','ymax', 'speed', 'FT', 'FToffset', ...
        'FToffsetRear', 'base_struct','sessions_included');

end

