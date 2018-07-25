function batch_align_pos(base_struct, reg_struct, varargin)
% batch_align_pos(base_struct, reg_struct, varargin)
%
% Aligns position data so that every session has the same bounds on the
% occupancy map and we can easily do correlations and other comparisons
% between firing across sessions.  Everything gets scaled to the
% trajectory/occupancy from the base session. IMPORTANT NOTE: running this
% WILL overwrite any previous Pos_align.mat files, so be careful (though
% not too careful since using this function is fairly quick). Note that
% this also aligns the position data to the imaging data found in
% FinalOutput.mat.
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
%       name_append: this will be appended to the Pos_align or
%       Pos_align_std_corr if you specify it.  Can also be a cell with
%       different names to append for each session
%
%       circ2square_use: logical with 1 indicating to use circle data
%       that has been transformed to square data (in Pos_trans.mat).  Will
%       save all data in Pos_align_trans.mat or
%       Pos_align_std_corr_trans.mat.
%
%   halfwindow: half the temporal smoothing window used - mainly obsolete,
%       but you may need to use if your data is offset for some reason.
%       Default = 0.
%
%   base_adjust: set to true the 1st time you run this or if you want to
%   overwrite a previously established alignment. True will overwrite and
%   pos_align.mat file found in the base session. default = false
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
%   AlignImagingToTracking
%    
%   'base_struct':
%   'sessions_included':
%
% NOTE: this may not work well for the 2 environment experiment since it
% does not account for any fish-eye distortions of the maze...should be
% good for most comparisons between the same mazes, however

close all;
global MasterDirectory;
MasterDirectory = 'C:\MasterData';

%% Parameters/default values
% Dump everything into one structure for future ease
base_struct = complete_struct(base_struct);
reg_struct = complete_struct(reg_struct);
sesh = [base_struct, reg_struct];
num_sessions = length(sesh);

p = inputParser;
p.addRequired('base_struct',@(x) isstruct(x));
p.addRequired('reg_struct',@(x) isstruct(x)); 
p.addParameter('manual_rot_overwrite',true,@(x) islogical(x));
p.addParameter('ratio_use',0.95,@(x) isscalar(x)); 
p.addParameter('auto_rotate_to_std',false,@(x) islogical(x) || x == 0 || x == 1); 
p.addParameter('manual_limits', false(1,num_sessions), @(x) islogical(x) ...
    && (length(x) == num_sessions || length(x) == 1));
p.addParameter('name_append','', @(x) ischar(x) || iscell(x) && length(x) == num_sessions); % Names to append to Pos_align file
p.addParameter('circ2square_use',false, @(x) islogical(x)); 
p.addParameter('TenaspisVer',4 ,@(x) isscalar(x) && x >= 3); 
p.addParameter('rotate_data', zeros(1,num_sessions), @(a) isnumeric(a) || ...
    islogical(a) && length(a) == num_sessions);
p.addParameter('skip_skew_fix', false, @islogical);
p.addParameter('suppress_output', false, @islogical);
p.addParameter('skip_trace_align', false, @islogical);
p.addParameter('base_adjust', false, @islogical); % Set to false if you want to incorporate new sessions into previously aligned ones - will leave base session alone.
p.addParameter('no_brain_data',false, @islogical);
p.addParameter('corner_align', false, @islogical); % option to NOT scale data but just align by the arena limits

p.parse(base_struct,reg_struct,varargin{:});
manual_rot_overwrite = p.Results.manual_rot_overwrite;
ratio_use = p.Results.ratio_use;
auto_rotate_to_std = p.Results.auto_rotate_to_std;
manual_limits = p.Results.manual_limits;
name_append = p.Results.name_append;
circ2square_use = p.Results.circ2square_use;
TenaspisVer = p.Results.TenaspisVer;
rotate_data = p.Results.rotate_data;
skip_skew_fix = p.Results.skip_skew_fix;
suppress_output = p.Results.suppress_output;
skip_trace_align = p.Results.skip_trace_align;
base_adjust = p.Results.base_adjust;
no_brain_data = p.Results.no_brain_data;
corner_align = p.Results.corner_align;
xmin = 10; ymin = 20; % Arbitrary point of reference - %%% NRK - need to adjust this for base_adjust = false - find95% and 95% limits of base data and adjust!!! 
xmin_ref = xmin; ymin_ref = ymin;

if length(manual_limits) == 1
    manual_limits = manual_limits*ones(1,num_sessions);
end
% Make name_append a cell
if ~iscell(name_append)
   temp = {name_append};
   clear name_append
   name_append = repmat(temp, 1, num_sessions);
end

dont_scale = false;
if corner_align
    dont_scale = true;
end

if base_adjust
    start_sesh = 1;
elseif ~base_adjust
    start_sesh = 2;
    base_align_file = fullfile(sesh(1).Location,['Pos_align' name_append{1} '.mat']);
end

curr_dir = cd;
%% 1: Load all sessions, and align to imaging data
pb = ProgressBar(num_sessions);
xcorner_all = nan(4,num_sessions); ycorner_all = nan(4,num_sessions);
for j = 1:num_sessions
    cd(sesh(j).Location);
    Pix2Cm = sesh(j).Pix2CM;
    if ~suppress_output
        disp(['Using ', num2str(Pix2Cm), ' as Pix2CM for ', sesh(j).Date, ' session ', num2str(sesh(j).Session)]);
    end
    
    if ~no_brain_data
    if TenaspisVer==4
        if ~suppress_output
            disp('Loading results from Tenaspis v4.');
        end
        HalfWindow = 0;
        load(fullfile(sesh(j).Location,'FinalOutput.mat'),'PSAbool','NeuronTraces');  
        LPtrace = NeuronTraces.LPtrace;
        DFDTtrace = NeuronTraces.DFDTtrace; 
        RawTrace = NeuronTraces.RawTrace;
        clear NeuronTraces;
    elseif TenaspisVer==3
        if ~suppress_output
            disp('Loading results from Tenaspis v3.');
        end
        HalfWindow = 0;
        load(fullfile(sesh(j).Location,'FinalOutput.mat'),'FT');
        load(fullfile(sesh(j).Location,'FinalTraces.mat'),'trace','difftrace','rawtrace');
        PSAbool = FT;
        LPtrace = trace;
        DFDTtrace = difftrace; 
        RawTrace = rawtrace;
        clear FT trace difftrace rawtrace;
    end
    
    % Align tracking and imaging
    [x,y,speed,PSAbool,FToffset,FToffsetRear,aviTime,time_interp,nframesinserted] = ...
        AlignImagingToTracking(Pix2Cm,PSAbool,HalfWindow, 'suppress_output', ...
        suppress_output);
    if ~skip_trace_align
        [~,~,~,LPtrace] = AlignImagingToTracking(Pix2Cm,LPtrace,HalfWindow, 'suppress_output', suppress_output);
        [~,~,~,DFDTtrace] = AlignImagingToTracking(Pix2Cm,DFDTtrace,HalfWindow, 'suppress_output', suppress_output);
        [~,~,~,RawTrace] = AlignImagingToTracking(Pix2Cm,RawTrace,HalfWindow, 'suppress_output', suppress_output);
    elseif skip_trace_align
        LPtrace = 'not performed';
        DFDTtrace = 'not performed';
        RawTrace = 'not performed';
    end
    elseif no_brain_data
        load(fullfile(sesh(j).Location,'Pos.mat'),'xpos_interp','ypos_interp',...
            'time_interp','MouseOnMazeFrame');  
        x = xpos_interp(MouseOnMazeFrame:end); %#ok<COLND>
        y = ypos_interp(MouseOnMazeFrame:end); %#ok<COLND>
        x = x.*Pix2Cm;
        y = y.*Pix2Cm;
        
        
        PSAbool = 'not performed';
        LPtrace = 'not performed';
        DFDTtrace = 'not performed';
        RawTrace = 'not performed';
        speed = 'not performed';
        FToffset = 'not performed';
        FToffsetRear = 'not performed';
        aviTime = 'not performed';
        nframesinserted = 'not performed';
    end
%%
% keyboard % NK - plot xcorner_all and ycorner_all on top of all data to check here!
%%
        
    % Set manual limits on data (e.g. to exclude wonky points outside the
    % maze)
    if ~manual_limits(j)
        x_for_limits = x;
        y_for_limits = y;
        sesh(j).ind_keep = true(1,length(x));

    elseif manual_limits(j)
        [x_for_limits, y_for_limits, sesh(j).ind_keep] = draw_manual_limits(...
            x, y, sesh(j).Env);
 
    end

    if corner_align
        [~, xAVI, yAVI, baseframe] = plot_traj2(sesh(j),'plot_vid', true,...
            'xy_append','AVI');
        BLok = false;
        while ~BLok
            title('Click on each corner of the arena')
            [xcor,ycor] = ginput(4); hold off
            imagesc(baseframe); set(gca,'YDir','Normal');
            hold on
            plot(xAVI,yAVI,'b-',xcor,ycor,'r*',mean(xcor),mean(ycor),'g*')
            okstr = input('is this ok (y/n)? : ', 's');
            BLok = strcmpi(okstr,'y');
        end
        % Need to scale xbl and ybl to x and y in cm, then subtract from x
        % and y.
%         avi_range = mean([range(xAVI(xAVI ~= 0)), range(yAVI(xAVI ~= 0))]);
        AVIgood = (xAVI > min(xcor) & xAVI < max(xcor) & ...
            yAVI > min(ycor) & yAVI < max(ycor));
        avi_range = mean([range(xAVI(AVIgood)), range(yAVI(AVIgood))]);
        cm_range = mean([range(x_for_limits(x_for_limits ~= 0)), ...
            range(y_for_limits(x_for_limits ~= 0))]);
        sf = cm_range/avi_range; % Should be equal to pix2CM...
        x_for_limits = x_for_limits - mean(xcor)*sf; 
        y_for_limits = y_for_limits - mean(ycor)*sf;
        xcorner_all(:,j) = (xcor - mean(xcor))*sf;
        ycorner_all(:,j) = (ycor - mean(ycor))*sf;
        
    end
    
   %% % Fix day-to-day mis-alignments in rotation of the maze
    if ~skip_skew_fix
        skewed = true;
        while skewed
            [rot_x,rot_y,rot_ang, rot_array, h] = rotate_traj(x_for_limits,...
                y_for_limits);
%             [rot_x,rot_y,rot_ang] = rotate_traj(x,y);
            plot(rot_x,rot_y);
            temp = (rot_array*[xcorner_all(:,j), ycorner_all(:,j)]')';
            
            hold on; plot(temp(:,1),temp(:,2),'r*')
            satisfied = input('Are you satisfied with the rotation? Enter y or n-->','s');
            skewed = ~strcmp(satisfied,'y');
            close(h)
        end
        xcorner_all(:,j) = temp(:,1);
        ycorner_all(:,j) = temp(:,2);
    elseif skip_skew_fix
        rot_x = x_for_limits;
        rot_y = y_for_limits;
%         rot_x = x;
%         rot_y = y;
        rot_ang = 'skipped';
        rot_array = eye(2);
    end
%%
%     % Auto-rotate back to standard configuration if indicated
%     if auto_rotate_to_std == 1
%         rot_corr = get_rot_from_db(sesh(j));
%         [x, y] = rotate_arena(x,y,rot_corr);
%     end
%     
%     if rotate_data(j) ~= 0
%         [x, y] = rotate_arena(x,y,rotate_data(j));
%     end
    
    sesh(j).x = x;
    sesh(j).y = y;
    sesh(j).PSAbool = PSAbool;
    sesh(j).LPtrace = LPtrace;
    sesh(j).DFDTtrace = DFDTtrace;
    sesh(j).RawTrace = RawTrace;
    sesh(j).speed = speed;
    sesh(j).FToffset = FToffset;
    sesh(j).FToffsetRear = FToffsetRear;
    
%     % Fix day-to-day mis-alignments in rotation of the maze
%     if ~skip_skew_fix
%         skewed = true;
%         while skewed
%             [rot_x,rot_y,rot_ang] = rotate_traj(x,y);
%             plot(rot_x,rot_y);
%             satisfied = input('Are you satisfied with the rotation? Enter y or n-->','s');
%             skewed = ~strcmp(satisfied,'y');
%         end
%     elseif skip_skew_fix
%         rot_x = x;
%         rot_y = y;
%         rot_ang = 'skipped';
%     end

    % Auto-rotate back to standard configuration if indicated
    if auto_rotate_to_std == 1
        rot_corr = get_rot_from_db(sesh(j));
        [rot_x, rot_y] = rotate_arena(x_for_limits, y_for_limits, rot_corr);
    end
    
    if rotate_data(j) ~= 0
        [rot_x, rot_y] = rotate_arena(x_for_limits, y_for_limits, rotate_data(j));
    end

    sesh(j).rot_x = rot_x;
    sesh(j).rot_y = rot_y;
    sesh(j).rot_ang = rot_ang;
    sesh(j).aviFrame = aviTime;
    sesh(j).time_interp = time_interp;
    sesh(j).nframesinserted = nframesinserted;
    
    pb.progress;
end

pb.stop;

cd(curr_dir)

%%
% keyboard

%% 2: Align position data for each session to the base session by using the 95% occupancy limits, save as Pos_align.mat
% Include base session in Pos_align for future reference

% Add in opportunity to manually select data limits to use here... if the
% flag you set for it is 1 - may want to do this by looking at the Notes
% section in MD


for j = 1:length(sesh)
if ~corner_align   
    x_for_limits = sesh(j).rot_x;
    y_for_limits = sesh(j).rot_y;
    %     if ~manual_limits(j)
    %         x_for_limits = sesh(j).rot_x;
    %         y_for_limits = sesh(j).rot_y;
    %         sesh(j).ind_keep = true(1,length(sesh(j).rot_x));
    %     elseif manual_limits(j)
    %         [x_for_limits, y_for_limits, sesh(j).ind_keep] = draw_manual_limits(...
    %             sesh(j).rot_x, sesh(j).rot_y, sesh(j).Env);
    %     end
    
    % Transform circle to square if indicated
    if circ2square_use && ~isempty(regexpi(sesh(j).Env,'octagon'))
        %Arena Size Parameters
        circle_radius = 14.33;
        square_side = 25.4;
        [x_for_limits, y_for_limits] = circ2square(x_for_limits, ...
            y_for_limits, square_side, circle_radius );
        sesh(j).rot_x = nan(size(sesh(j).rot_x));
        sesh(j).rot_y = nan(size(sesh(j).rot_y));
        sesh(j).rot_x(sesh(j).ind_keep) = x_for_limits;
        sesh(j).rot_y(sesh(j).ind_keep) = y_for_limits;
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
    
    % Re-do the above but using raw x and y values adjusted by corners only.
elseif corner_align
    sesh(j).x_adj = sesh(j).rot_x + xmin;
    sesh(j).y_adj = sesh(j).rot_y + ymin;
end

% Adjust xmin and ymin if adding in new data to existing aligned
% sessions(s).
if j == 1 && ~base_adjust %% NRK - check here!
    load(base_align_file,'x_adj_cm', 'y_adj_cm');
    xmin = xmin + min(x_adj_cm) - min(sesh(1).x_adj);
    ymin = ymin + min(y_adj_cm) - min(sesh(1).y_adj);
    
    sesh(1).x_adj = x_adj_cm;
    sesh(1).y_adj = y_adj_cm;
end

end

% keyboard

%% 4: Concatenate ALL position data into one X and one Y vector, and get Xedges and Yedges based on this
x_all = [];
y_all = [];
for j = 1:length(sesh)
    x_all = [x_all sesh(j).x_adj]; % Only use data within the limits you drew for this!
    y_all = [y_all sesh(j).y_adj]; % Only use data within the limits you drew for this!
%     x_all = [x_all sesh(j).x_adj(sesh(j).ind_keep)]; % Only use data within the limits you drew for this!
%     y_all = [y_all sesh(j).y_adj(sesh(j).ind_keep)]; % Only use data within the limits you drew for this!
end

%% 5: Get xmin, xmax, ymin, and ymax

if base_adjust
    if any(isnan(xcorner_all(:)))
        xmax = max(x_all)+1; % give yourself a buffer of 1 cm on each side just in case you want to align future sessions
        xmin = min(x_all)-1;
        ymax = max(y_all)+1;
        ymin = min(y_all)-1;
    elseif all(~isnan(xcorner_all(:)))
        xmax = max(xcorner_all(:))+1 + xmin_ref; % give yourself a buffer of 1 cm on each side just in case you want to align future sessions
        xmin = min(xcorner_all(:))-1 + xmin_ref;
        ymax = max(ycorner_all(:))+1 + ymin_ref;
        ymin = min(ycorner_all(:))-1 + ymin_ref;
    end
elseif ~base_adjust % Set limits to those from base session if not adjusting the base session.
    load(base_align_file,'xmin', 'xmax', 'ymin', 'ymax');
end

% keyboard

%% 6: Save Xedges, Yedges in base session for future reference along with all sessions aligned to it.
% Also save adjusted position data for future use...

sessions_included = [base_struct reg_struct];

for j = start_sesh:length(sesh)
    if ~auto_rotate_to_std%     x_adj_cm = sesh(j).x_adj;
%     y_adj_cm = sesh(j).y_adj;
%     speed = sesh(j).speed;
%     PSAbool = sesh(j).PSAbool;
%     LPtrace = sesh(j).LPtrace;
%     DFDTtrace = sesh(j).DFDTtrace;
%     RawTrace = sesh(j).RawTrace;
%     FToffset = sesh(j).FToffset;
%     FToffsetRear = sesh(j).FToffsetRear;
%     aviFrame = sesh(j).aviFrame;
%     time_interp = sesh(j).time_interp;
%     nframesinserted = sesh(j).nframesinserted;
%     Pix2CM = sesh(j).Pix2CM;
    ik = sesh(j).ind_keep;
    x_adj_cm = sesh(j).x_adj;
    y_adj_cm = sesh(j).y_adj;
    speed = sesh(j).speed(ik);
    PSAbool = sesh(j).PSAbool(:,ik);
    LPtrace = sesh(j).LPtrace(:,ik);
    DFDTtrace = sesh(j).DFDTtrace(:,ik);
    RawTrace = sesh(j).RawTrace(:,ik);
    FToffset = sesh(j).FToffset;
    FToffsetRear = sesh(j).FToffsetRear;
    aviFrame = sesh(j).aviFrame(ik);
    time_interp = sesh(j).time_interp(ik);
    nframesinserted = sesh(j).nframesinserted;
    Pix2CM = sesh(j).Pix2CM;
    if find(~ik,1,'first') <= FToffset % Un-checked bugfix
        FToffset = FToffset - length(find(~ik) <= FToffset);
        disp('chopped frames before FToffset - debug here')
        keyboard
    end

    save(fullfile(sesh(j).Location,['Pos_align' name_append{j} '.mat']),...
        'x_adj_cm','y_adj_cm','xmin','xmax','ymin','ymax','speed',...
        'PSAbool','LPtrace','DFDTtrace','RawTrace','FToffset',...
        'nframesinserted','time_interp','FToffsetRear','aviFrame',...
        'base_struct','sessions_included','auto_rotate_to_std', 'Pix2CM');
    elseif auto_rotate_to_std
        % finish here - save as a different filename?
        save(fullfile(sesh(j).Location,...
            ['Pos_align_std_corr' name_append{j} '.mat']),'x_adj_cm',...
            'y_adj_cm','xmin','xmax','ymin','ymax','speed','PSAbool',...
            'LPtrace','DFDTtrace','RawTrace','FToffset','nframesinserted',...
            'time_interp','FToffsetRear','aviFrame','base_struct',...
            'sessions_included', 'auto_rotate_to_std', 'Pix2CM');
    end
end

%% 7: Plot everything as a check
figure(100)
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

%% Fill-in partially completely base_struct or reg_struct
function [struct_out] = complete_struct(struct_in)
num_sessions = length(struct_in);

for j = 1:num_sessions
    [~, struct_out(j)] = ChangeDirectory(struct_in(j).Animal,...
        struct_in(j).Date, struct_in(j).Session,0);
end

end
