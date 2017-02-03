function [xpos_interp,ypos_interp,start_time,MoMtime,time_interp,AVItime_interp,nframesinserted] = PreProcessMousePosition_auto(filepath, auto_thresh,varargin)
% [xpos_interp,ypos_interp,start_time,MoMtime] = PreProcessMousePosition_auto(filepath, auto_thresh,...)
% Function to correct errors in mouse tracking.  Runs once through the
% entire sessions automatically having you edit any events above a velocity
% threshold (set by 'auto_thresh', suggest setting this to 0.01 or 0.02).
%
% INPUTS
%   filepath: pathname to DVT file. Must reside in the same directory as
%   the AVI file it matches, and there must be only ONE DVT and ONE AVI
%   file in this directory
% 
%   auto_thresh: proportion of timestamps you wish to edit - 0.01 will have
%   you edit all timestamps where the velocity of the mouse is in the top
%   1% of the distribution of all velocities - suggest starting at 0.01 or
%   0.02 for a typical session, but stepping down to 0.001 or smaller if
%   you do a second pass.
%
%   'update_pos_realtime' (optional): Default is 0. Set to 1 if you want to watch the
%   position getting updated with each click of the mouse, but I don't
%   suggest it because it tends to cause weird crashes when MATLAB can't
%   figure out which figure it should actually be plotting to.
%
%   'epoch_length_lim': will not auto-correct any epochs over this length
%   where the mouse is at 0,0 or above the velocity threhold - suggest
%   using if the mouse is off the maze for a long time.
%
%   OUTPUTS (all saved in Pos.mat, along with some others)
%   xpos_interp, ypos_interp: smoothed, corrected position data
%   interpolated to match the frame rate of the imaging data (hardcoded at
%   20 fps)
%
%   start_time: start of DVT file
%
%   MoMtime: the time that the mouse starts running on the maze

close all;

%% Get varargin

update_pos_realtime = 0; % Default setting
epoch_length_lim = 200; % default
for j = 1:length(varargin)
   if strcmpi('update_pos_realtime', varargin{j})
      update_pos_realtime = varargin{j+1};
   end
   if strcmpi('epoch_length_lim', varargin{j})
      epoch_length_lim = varargin{j+1};
   end
   
end

%%
% Script to take position data at given timestamps and output and interpolate 
% to any given timestamps.

PosSR = 30; % native sampling rate in Hz of position data (used only in smoothing)
aviSR = 30.0003; % the framerate that the .avi thinks it's at
cluster_thresh = 40; % For auto thresholding - any time there are events above
% the velocity threshold specified by auto_thresh that are less than this
% number of frames apart they will be grouped together


% Import position data from DVT file
try
pos_data = importdata(filepath);
%f_s = max(regexp(filepath,'\'))+1;
%mouse_name = filepath(f_s:f_s+2);
%date = [filepath(f_s+3:f_s+4) '-' filepath(f_s+5:f_s+6) '-' filepath(f_s+7:f_s+10)];

% Parse out into invididual variables
frame = pos_data(:,1);
time = pos_data(:,2); % time in seconds
Xpix = pos_data(:,3); % x position in pixels (can be adjusted to cm)
Ypix = pos_data(:,4); % y position in pixels (can be adjusted to cm)
catch
% Video.txt is there instead of Video.DVT
pos_data = importdata('Video.txt');
Xpix = pos_data.data(:,6);
Ypix = pos_data.data(:,7);
time = pos_data.data(:,4);
end

%WM Edit: Check for correct Cineplex sampling rate 
dt = [0.03; round(diff(time),2)]; 
bad = dt > 0.06; % Look for anything close to or greater than 2x the sampling rate 
bad_time = time;
cum_l = 0;
if any(bad)
    disp('Dropped frames detected! Correcting as best I can...'); 
    im = regionprops(bad,'Area','PixelIdxList'); 
    
    for i=1:length(im)
        if im(i).Area > 1
            s = im(i).PixelIdxList(1)+cum_l;
            e = im(i).PixelIdxList(end)+cum_l;

            if s < length(time)
                tInsert = [(time(s)+1/aviSR):(1/aviSR):(time(e)-1/aviSR)]';
                l = length(tInsert);  
                cum_l = cum_l + l;
                XpixInsert = Xpix(s-1)*ones(l,1); 
                YpixInsert = Ypix(s-1)*ones(l,1);

                time = [time(1:s); tInsert; time(e:end)];
                Xpix = [Xpix(1:s); XpixInsert; Xpix(e:end)];
                Ypix = [Ypix(1:s); YpixInsert; Ypix(e:end)];
            end
        end       
    end
end


xAVI = Xpix*.6246;
yAVI = Ypix*.6246;

figure(777);plot(Xpix,Ypix);title('pre-corrected data');

try
    h1 = implay('Raw.AVI');
    obj = VideoReader('Raw.AVI');
catch
    avi_filepath = ls('*.avi');
    h1 = implay(avi_filepath);
    disp(['Using ' avi_filepath ])
    obj = VideoReader(avi_filepath);
end

if exist('Pos_temp.mat','file') || exist('Pos.mat','file')
    % Determine if either Pos_temp or Pos file already exists in the
    % directory, and prompt user to load it up if they want to continue
    % editing it.
    if exist('Pos_temp.mat','file') && ~exist('Pos.mat','file')
        use_temp = input('Pos_temp.mat detected.  Enter "y" to use or "n" to start from scratch: ','s');
        load_file = 'Pos_temp.mat';
    elseif exist('Pos.mat','file')
        use_temp = input('Previous Pos.mat detected.  Enter "y" to use or "n" to start from scratch: ','s');
        load_file = 'Pos.mat';
    end
    if strcmpi(use_temp,'y')
        load(load_file,'Xpix', 'Ypix', 'xAVI', 'yAVI', 'MoMtime', 'MouseOnMazeFrame');
        MoMtime
    else
        MouseOnMazeFrame = input('on what frame number does Mr. Mouse arrive on the maze??? --->');
        MoMtime = MouseOnMazeFrame/aviSR+time(1)
    end
else
    MouseOnMazeFrame = input('on what frame number does Mr. Mouse arrive on the maze??? --->');
    MoMtime = MouseOnMazeFrame/aviSR+time(1)
end
close(h1); % Close Video Player

% Get initial velocity profile for auto-thresholding
vel_init = hypot(diff(Xpix),diff(Ypix))/(time(2)-time(1));
vel_init = [vel_init; vel_init(end)];
% vel_init = [vel_init(1); vel_init];
[fv, xv] = ecdf(vel_init);
if exist('auto_thresh','var')
    auto_vel_thresh = min(xv(fv > (1-auto_thresh)));
else
    auto_vel_thresh = max(vel_init)+1;
    auto_thresh = nan; % Don't perform any autocorrection if not specified
end

% start auto-correction of anything above threshold
auto_frames = (Xpix == 0 | Ypix == 0 | vel_init > auto_vel_thresh) & time > MoMtime;

% Determine if auto thresholding applies
if sum(auto_frames) > 0 && ~isnan(auto_thresh)
    auto_thresh_flag = 1;
    [ on, off ] = get_on_off( auto_frames );
    [ epoch_start, epoch_end ] = cluster_epochs( on, off, cluster_thresh );
    n_epochs = length(epoch_start);
    
    % Apply epoch length limits if applicable
    if ~isempty(epoch_length_lim)
       epoch_lengths = epoch_end - epoch_start;
       epoch_start = epoch_start(epoch_lengths < epoch_length_lim);
       epoch_end = epoch_end(epoch_lengths < epoch_length_lim);
       n_epochs = length(epoch_start);
    end
else %if sum(auto_frames) == 0
    auto_thresh_flag = 0;
end

% keyboard

figure(555);
hx0 = subplot(4,3,1:3);plot(time,Xpix);xlabel('time (sec)');ylabel('x position (cm)');yl = get(gca,'YLim');line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');axis tight;
hy0 = subplot(4,3,4:6);plot(time,Ypix);xlabel('time (sec)');ylabel('y position (cm)');yl = get(gca,'YLim');line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');axis tight;
linkaxes([hx0 hy0],'x');

v0 = readFrame(obj);
MorePoints = 'y';
%length(time);

%Draw a mask for the maze.
figure;
imagesc(flipud(v0)); title('Draw a mask for the maze'); 
maze = roipoly; 

n = 1; %first_time = 1;
while (strcmp(MorePoints,'y')) || strcmp(MorePoints,'m') || isempty(MorePoints)
%     if first_time == 1
%         hx0 = subplot(4,3,1:3); plot(time,Xpix); xlabel('time (sec)'); ylabel('x position (cm)');
%         hold on;yl = get(gca,'YLim');line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');hold off;axis tight;
%         hy0 = subplot(4,3,4:6); plot(time,Ypix); xlabel('time (sec)'); ylabel('y position (cm)');
%         hold on;yl = get(gca,'YLim');line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');hold off;axis tight;
%         first_time = 0;
%         linkaxes([hx0 hy0],'x');
%     end
  if auto_thresh_flag == 0 || isempty(epoch_start)
      figure(555);
      MorePoints = input('Is there a flaw that needs to be corrected?  [y/n/manual correct (m)] -->','s');
  else
      MorePoints = 'y'; pause(1)
  end

  
  if strcmp(MorePoints,'y')
      if auto_thresh_flag == 0 || isempty(epoch_start)
          FrameSelOK = 0;
          while (FrameSelOK == 0)
              display('click on the good points around the flaw then hit enter');
              [DVTsec,~] = ginput(2); % DVTsec is start and end time in DVT seconds
              sFrame = findclosest(time,DVTsec(1)); % index of start frame
              eFrame = findclosest(time,DVTsec(2)); % index of end frame
              aviSR*sFrame;
              
              if (sFrame/aviSR > obj.Duration || eFrame/aviSR > obj.Duration)
                  
                  continue;
              end
%               obj.currentTime = sFrame/aviSR; % sFrame is the correct frame #, but .avi reads are done according to time
%               v = readFrame(obj);
              FrameSelOK = 1;
              
          end
          
      elseif auto_thresh_flag == 1 % Input times from auto_threholded vector
          sFrame = max([1 epoch_start(n)- 6]);
          eFrame = min([length(time) epoch_end(n) + 6]);
          
          % Turn on manual thresholding once you correct all epochs above
          % the velocity threshold
          if n == n_epochs
              auto_thresh_flag = 0;
          else
              n = n + 1;
          end
      end
    obj.currentTime = sFrame/aviSR; % sFrame is the correct frame #, but .avi reads are done according to time
    v = readFrame(obj);
    
    framesToCorrect = sFrame:eFrame;
    if eFrame >= max(time)
        framesToCorrect = sFrame:eFrame-2; % Fix case where last frame needs to be corrected
    end
    frame_use_index = 1:floor(length(framesToCorrect)/2);
    frame_use_num = length(frame_use_index);
    
    edit_start_time = time(sFrame);
    edit_end_time = time(eFrame);
    
    % Set marker colors to be green for the first 1/3, yellow for the 2nd
    % 1/3, and red for the final 1/3
    marker = {'go' 'yo' 'ro'};
    marker_face = {'g' 'y' 'r'};
    marker_fr = ones(size(frame_use_index));
    num_markers = size(marker,2);
    for jj = 1:num_markers-1
        marker_fr(floor(jj*frame_use_num/num_markers)+1:...
            floor((jj+1)*frame_use_num/num_markers)) = ...
            (jj+1)*ones(size(floor(jj*frame_use_num/num_markers)+1:...
            floor((jj+1)*frame_use_num/num_markers)));
    end
   
    
    disp(['You are currently editing from ' num2str(edit_start_time) ...
        ' sec to ' num2str(edit_end_time) ' sec.'])
     
    for i = frame_use_index
        
        if update_pos_realtime == 1
            figure(555)
            % Plot updated coordinates and velocity
            % plot the current sub-trajectory
            subplot(4,3,11);
            imagesc(flipud(v));hold on;
            plot(xAVI(sFrame:eFrame),yAVI(sFrame:eFrame),'LineWidth',1.5);hold off;title('chosen segment');
            
            % plot the current total trajectory
            subplot(4,3,10);
            imagesc(flipud(v));hold on;
            plot(xAVI(MouseOnMazeFrame:end),yAVI(MouseOnMazeFrame:end),'LineWidth',1.5);hold off;title('overall trajectory (post mouse arrival)');
        end
        
        % plot the current video frame
        %framesToCorrect(i*2);
        obj.currentTime = framesToCorrect(i*2)/aviSR;
        v = readFrame(obj);
        figure(1702);pause(0.01);
        gcf;
        imagesc(flipud(v));title('click here');
        % plot the existing position marker on top
        hold on;plot(xAVI(sFrame+i*2),yAVI(sFrame+i*2),marker{marker_fr(i)},'MarkerSize',4);
%         display(['Time is ' num2str(time(sFrame+i*2)) ' seconds. Click the mouse''s back']);
        
        %Subtract current frame from reference, then flip and smooth. Next,
        %run regionprops. 
        d = imgaussfilt(flipud(rgb2gray(v0-v)),10); 
        stats = regionprops(d>18 & maze,'area','centroid','majoraxislength','minoraxislength');        
        
        %Find the blob that corresponds to the mouse. 
        MouseBlob = find(   [stats.Area] > 200 & ...
                            [stats.MajorAxisLength] > 10 & ...
                            [stats.MinorAxisLength] > 10);
        if length(MouseBlob)==1     
            xm = stats(MouseBlob).Centroid(1); 
            ym = stats(MouseBlob).Centroid(2);
        elseif length(MouseBlob)>1
            %Get mouse position on the previous frame. 
            previousX = xAVI(sFrame+i*2-1);
            previousY = yAVI(sFrame+i*2-1); 
            
            %Possible mouse blobs. 
            putativeMouse = [stats(MouseBlob).Centroid];
            putativeMouseX = putativeMouse(1:2:end);
            putativeMouseY = putativeMouse(2:2:end); 
            
            %Find which blob is closest to the mouse's location in the
            %previous frame. 
            whichMouseX = findclosest(putativeMouseX,previousX); 
            whichMouseY = findclosest(putativeMouseY,previousY); 
            
            %If they agree, use that blob. 
            if whichMouseX == whichMouseY
                xm = stats(MouseBlob(whichMouseX)).Centroid(1);
                ym = stats(MouseBlob(whichMouseY)).Centroid(2); 
            else
                %keyboard;
                [xm,ym] = ginput(1); 
            end
        else          
            %keyboard;
            [xm,ym] = ginput(1);
        end
        
        % apply corrected position to current frame
        xAVI(sFrame+i*2) = xm;
        yAVI(sFrame+i*2) = ym;
        Xpix(sFrame+i*2) = ceil(xm/0.6246);
        Ypix(sFrame+i*2) = ceil(ym/0.6246);
        
        % interpolate and apply correct position for previous frame
        xAVI(sFrame+i*2-1) = xAVI(sFrame+i*2-2)+(xm-xAVI(sFrame+i*2-2))/2;
        yAVI(sFrame+i*2-1) = yAVI(sFrame+i*2-2)+(ym-yAVI(sFrame+i*2-2))/2;
        Xpix(sFrame+i*2-1) = ceil(xAVI(sFrame+i*2-1)/0.6246);
        Ypix(sFrame+i*2-1) = ceil(yAVI(sFrame+i*2-1)/0.6246);
        
        
        % plot marker
        plot(xm,ym,marker{marker_fr(i)},'MarkerSize',4,'MarkerFaceColor',marker_face{marker_fr(i)});hold off;
    end
    disp(['You just edited from ' num2str(edit_start_time) ...
        ' sec to ' num2str(edit_end_time) ' sec.']);
   
    close(1702);
    
    % plot updated velocity
    figure(555);
    subplot(4,3,7:9);
    vel = hypot(diff(Xpix),diff(Ypix))/(time(2)-time(1));
    vel = [vel; vel(end)]; % Make the vectors the same size
    plot(time(MouseOnMazeFrame:end),vel(MouseOnMazeFrame:end));
    hold on
    plot(time([sFrame eFrame]),vel([sFrame eFrame]),'ro'); % plot start and end points of last edit
    if auto_thresh_flag == 1
        % Get indices for all remaining times that fall above the auto 
        % threshold that have not been corrected
        ind_red = auto_frames & time > time(eFrame); 
        hold on
        plot(time(ind_red),vel(ind_red),'ro');
        hold off
    end
    hold off;axis tight;xlabel('time (sec)');ylabel('velocity (units/sec)'); 
    set(gca,'xtick',[]); 
    xlim_use = get(gca,'XLim'); hv = gca;
    
    % plot updated x and y values
    hx = subplot(4,3,1:3); plot(time,Xpix); hold on; 
    plot(time([sFrame eFrame]),Xpix([sFrame eFrame]),'ro'); % plot start and end points of last edit
    xlabel('time (sec)'); ylabel('x position (cm)');
    hold on; yl = get(gca,'YLim'); line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');
    hold off; axis tight; % set(gca,'XLim',[sFrame/aviSR eFrame/aviSR]); hx = gca;
    
    hy = subplot(4,3,4:6); plot(time,Ypix); hold on; 
    plot(time([sFrame eFrame]),Ypix([sFrame eFrame]),'ro'); % plot start and end points of last edit
    xlabel('time (sec)'); ylabel('y position (cm)');
    hold on; yl = get(gca,'YLim'); line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');
    hold off; axis tight; % set(gca,'XLim',[sFrame/aviSR eFrame/aviSR]); hy = gca;
    
    linkaxes([hx, hy, hv],'x'); % Link axes zooms along time dimension together
    
    drawnow % Make sure everything gets updated properly!
    
    % NRK edit
    save Pos_temp.mat Xpix Ypix xAVI yAVI MoMtime MouseOnMazeFrame
    
  continue
  end
  
  if (strcmp(MorePoints,'g'))
      % generate a movie and show it
      for i = 1:length(time)
        obj.currentTime = i/aviSR; % sFrame is the correct frame #, but .avi reads are done according to time
        v = readFrame(obj);
        figure(6156);
        imagesc(flipud(v));hold on;
        plot(xAVI(i),yAVI(i),'or','MarkerSize',5,'MarkerFaceColor','r');hold off;
        F(i) = getframe(gcf);
      end
      save F.mat F;implay(F);pause;
  end
  
if strcmp(MorePoints,'m')
      %Copied from above. Frame selection. 
      FrameSelOK = 0;
          while (FrameSelOK == 0)
              display('click on the good points around the flaw then hit enter');
              [DVTsec,~] = ginput(2); % DVTsec is start and end time in DVT seconds
              sFrame = findclosest(time,DVTsec(1)); % index of start frame
              eFrame = findclosest(time,DVTsec(2)); % index of end frame
              aviSR*sFrame;
              
              if (sFrame/aviSR > obj.Duration || eFrame/aviSR > obj.Duration)
                  
                  continue;
              end
%               obj.currentTime = sFrame/aviSR; % sFrame is the correct frame #, but .avi reads are done according to time
%               v = readFrame(obj);
              FrameSelOK = 1;
              
          end
          
    %Copied from above. Old way of correcting video.     
    obj.currentTime = sFrame/aviSR; % sFrame is the correct frame #, but .avi reads are done according to time
    v = readFrame(obj);
       
    framesToCorrect = sFrame:eFrame;
    if eFrame >= max(time)
        framesToCorrect = sFrame:eFrame-2; % Fix case where last frame needs to be corrected
    end
    frame_use_index = 1:floor(length(framesToCorrect)/2);
    frame_use_num = length(frame_use_index);
    
    edit_start_time = time(sFrame);
    edit_end_time = time(eFrame);
    
    % Set marker colors to be green for the first 1/3, yellow for the 2nd
    % 1/3, and red for the final 1/3
    marker = {'go' 'yo' 'ro'};
    marker_face = {'g' 'y' 'r'};
    marker_fr = ones(size(frame_use_index));
    num_markers = size(marker,2);
    for jj = 1:num_markers-1
        marker_fr(floor(jj*frame_use_num/num_markers)+1:...
            floor((jj+1)*frame_use_num/num_markers)) = ...
            (jj+1)*ones(size(floor(jj*frame_use_num/num_markers)+1:...
            floor((jj+1)*frame_use_num/num_markers)));
    end
   
    
    disp(['You are currently editing from ' num2str(edit_start_time) ...
        ' sec to ' num2str(edit_end_time) ' sec.'])
     
    for i = frame_use_index
        
        if update_pos_realtime == 1
            figure(555)
            % Plot updated coordinates and velocity
            % plot the current sub-trajectory
            subplot(4,3,11);
            imagesc(flipud(v));hold on;
            plot(xAVI(sFrame:eFrame),yAVI(sFrame:eFrame),'LineWidth',1.5);hold off;title('chosen segment');
            
            % plot the current total trajectory
            subplot(4,3,10);
            imagesc(flipud(v));hold on;
            plot(xAVI(MouseOnMazeFrame:end),yAVI(MouseOnMazeFrame:end),'LineWidth',1.5);hold off;title('overall trajectory (post mouse arrival)');
        end
        
        % plot the current video frame
        framesToCorrect(i*2);
        obj.currentTime = framesToCorrect(i*2)/aviSR;
        v = readFrame(obj);
        figure(1702);pause(0.01);
        gcf;
        imagesc(flipud(v));title('click here');
        % plot the existing position marker on top
        hold on;plot(xAVI(sFrame+i*2),yAVI(sFrame+i*2),marker{marker_fr(i)},'MarkerSize',4);
     
        %Correct frames here!
        [xm,ym] = ginput(1);
        
        % apply corrected position to current frame
        xAVI(sFrame+i*2) = xm;
        yAVI(sFrame+i*2) = ym;
        Xpix(sFrame+i*2) = ceil(xm/0.6246);
        Ypix(sFrame+i*2) = ceil(ym/0.6246);
        
        % interpolate and apply correct position for previous frame
        xAVI(sFrame+i*2-1) = xAVI(sFrame+i*2-2)+(xm-xAVI(sFrame+i*2-2))/2;
        yAVI(sFrame+i*2-1) = yAVI(sFrame+i*2-2)+(ym-yAVI(sFrame+i*2-2))/2;
        Xpix(sFrame+i*2-1) = ceil(xAVI(sFrame+i*2-1)/0.6246);
        Ypix(sFrame+i*2-1) = ceil(yAVI(sFrame+i*2-1)/0.6246);
        
        
        % plot marker
        plot(xm,ym,marker{marker_fr(i)},'MarkerSize',4,'MarkerFaceColor',marker_face{marker_fr(i)});hold off;
    end
    disp(['You just edited from ' num2str(edit_start_time) ...
        ' sec to ' num2str(edit_end_time) ' sec.']);
   
    close(1702);
    
    % plot updated velocity
    figure(555);
    subplot(4,3,7:9);
    vel = sqrt(diff(Xpix).^2+diff(Ypix).^2)/(time(2)-time(1));
    vel = [vel; vel(end)]; % Make the vectors the same size
    plot(time(MouseOnMazeFrame:end),vel(MouseOnMazeFrame:end));
    hold on
    plot(time([sFrame eFrame]),vel([sFrame eFrame]),'ro'); % plot start and end points of last edit
    if auto_thresh_flag == 1
        % Get indices for all remaining times that fall above the auto 
        % threshold that have not been corrected
        ind_red = auto_frames & time > time(eFrame); 
        hold on
        plot(time(ind_red),vel(ind_red),'ro');
        hold off
    end
    hold off;axis tight;xlabel('time (sec)');ylabel('velocity (units/sec)');
    xlim_use = get(gca,'XLim'); hv = gca;
    
    % plot updated x and y values
    hx = subplot(4,3,1:3); plot(time,Xpix); hold on; 
    plot(time([sFrame eFrame]),Xpix([sFrame eFrame]),'ro'); % plot start and end points of last edit
    xlabel('time (sec)'); ylabel('x position (cm)');
    hold on; yl = get(gca,'YLim'); line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');
    hold off; axis tight; % set(gca,'XLim',[sFrame/aviSR eFrame/aviSR]); hx = gca;
    
    hy = subplot(4,3,4:6); plot(time,Ypix); hold on; 
    plot(time([sFrame eFrame]),Ypix([sFrame eFrame]),'ro'); % plot start and end points of last edit
    xlabel('time (sec)'); ylabel('y position (cm)');
    hold on; yl = get(gca,'YLim'); line([MoMtime MoMtime], [yl(1) yl(2)],'Color','r');
    hold off; axis tight; % set(gca,'XLim',[sFrame/aviSR eFrame/aviSR]); hy = gca;
    
    linkaxes([hx, hy, hv],'x'); % Link axes zooms along time dimension together
    
    drawnow % Make sure everything gets updated properly!
    
    % NRK edit
    save Pos_temp.mat Xpix Ypix xAVI yAVI MoMtime MouseOnMazeFrame
    
  continue
end
end

Xpix_filt = NP_QuickFilt(Xpix,0.0000001,1,PosSR);
Ypix_filt = NP_QuickFilt(Ypix,0.0000001,1,PosSR);

% if size(pos_data,2) == 5
%     motion = pos_data(:,5);
% end

AVIobjTime = zeros(1,length(time)); 
for i = 1:length(time)
    AVIobjTime(i) = i/aviSR;
end

frame_rate_emp = round(1/mean(diff(time))); % empirical frame rate (frames/sec)

% Generate times to match brain imaging data timestamps
fps_brainimage = 20; % frames/sec for brain image timestamps

start_time = ceil(min(time)*fps_brainimage)/fps_brainimage;
max_time = floor(max(time)*fps_brainimage)/fps_brainimage;
time_interp = start_time:1/fps_brainimage:max_time;                         %20 Hz timer starts when you hit record.

if (max(time_interp) >= max_time)
    time_interp = time_interp(1:end-1);
end

%% Do Linear Interpolation

% Get appropriate time points to interpolate for each timestamp
time_index = arrayfun(@(a) [max(find(a >= time)) min(find(a < time))],...
    time_interp,'UniformOutput',0);
time_test_cell = arrayfun(@(a) a,time_interp,'UniformOutput',0);

xpos_interp = cellfun(@(a,b) lin_interp(time(a), Xpix_filt(a),...           %20 Hz timer starts when you hit record.
    b),time_index,time_test_cell);

ypos_interp = cellfun(@(a,b) lin_interp(time(a), Ypix_filt(a),...           %20 Hz timer starts when you hit record.
    b),time_index,time_test_cell);  

AVItime_interp = cellfun(@(a,b) lin_interp(time(a), AVIobjTime(a),...       %20 Hz timer starts when you hit record.
    b),time_index,time_test_cell);

nframesinserted = round(length(AVItime_interp) - length(bad_time)*0.6664);

% Save all filtered data as well as raw data in case you want to go back
% and fix an error you discover later on
save Pos.mat xpos_interp ypos_interp time_interp start_time MoMtime Xpix Ypix xAVI yAVI MouseOnMazeFrame AVItime_interp nframesinserted;

end