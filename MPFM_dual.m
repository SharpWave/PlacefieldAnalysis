function [ ] = MPFM_dual(out_avifile,infile,clims,Pix2Cm,varargin)
% MPFM_dual(out_avifile,infile,clims,Pix2Cm)
%   make a fun movie of our mouse and his placefields
%
% INPUTS
%   out_avifile: name of the output AVI file you wish to create
%
%   infile: h5 brain imaging file you wish to show
%
%   clims: ?
%
%   Pix2Cm: conversion factor, should be whatever you used for running
%   CalculatePlacefields
%
%   varargins:
%       'brain_only': plot brain imaging only
%       'PF_only': plot placefields on video tracking only

video_type = 1; % default to show both videos side-by-side
if length(varargin) == 1
   if strcmpi(varargin{1},'brain_only')
       video_type = 2;
   elseif strcmpi(varargin{1},'PF_only')
       video_type = 3;
   end
end

close all;

if (~exist('Pix2Cm'))
    Pix2Cm = 0.15;
    display('assuming room 201b');
    % factor for 201a is 0.0709
    % else
    %     if (strcmp(RoomStr,'201a'))
    %         Pix2Cm = 0.0709;
    %         display('Room 201a');
    %     end
end

aviSR = 30.0003;

try
    %h1 = implay('Raw.AVI');
    obj = VideoReader('Raw.AVI');
catch
    avi_filepath = ls('*.avi');
    %h1 = implay(avi_filepath);
    disp(['Using ' avi_filepath ])
    obj = VideoReader(avi_filepath);
end

load PlaceMaps.mat;
load PFstats.mat;


NumFrames = length(x);
NumNeurons = length(NeuronImage);
Xdim = size(NeuronImage{1},1);
Ydim = size(NeuronImage{1},2);

figure;
set(gcf,'Position',[534 72 1171 921]);

aviobj = VideoWriter(out_avifile,'MPEG-4');
aviobj.FrameRate = 20;
aviobj.Quality = 90;
open(aviobj);

% assign each neuron a color
colors = rand(NumNeurons,3);

% convert x and y to .avi coords
xAVI = x/Pix2Cm*0.625;
yAVI = y/Pix2Cm*0.625;

% convert Xbin and Ybin to x and y
Xd = Xedges(2)-Xedges(1);
Yd = Yedges(2)-Yedges(1);
if video_type == 1 || video_type == 3

end

for i = 1:length(Xedges)
    Xb2AVI(i) = (Xedges(i)+Xd/2)/Pix2Cm*0.625;
end

for i = 1:length(Yedges)
    Yb2AVI(i) = (Yedges(i)+Yd/2)/Pix2Cm*0.625;
end

% for each neuron
try
for j = 1:NumNeurons
    % get PF outline (if avail)
    WhichField = MaxPF(j);
    temp = zeros(size(TMap{1}));
    tp = PFpixels{j,WhichField};
    temp(tp) = 1;
    nt(j) = size(NP_FindSupraThresholdEpochs(FT(j,:),eps),1);
    % plot PF outline (using correct color)
    b = bwboundaries(temp,4);
    
    if(~isempty(b))
        yt{j} = Yb2AVI(b{1}(:,2));
        xt{j} = Xb2AVI(b{1}(:,1));
        xt{j}= xt{j}+(rand(size(xt{j}))-0.5)/2;
        yt{j}= yt{j}+(rand(size(yt{j}))-0.5)/2;
        %colors(j,:)
        %plot(xt,yt,'Color',colors(j,:),'LineWidth',5);
    end
end
catch
    disp('MPFM_dual error catching')
    keyboard
end



for i = 1:NumFrames
    
    % load correct Plexon movie frame
    % calculate correct frame based on iteration and offsets
    obj.currentTime = aviFrame(i);
    v = readFrame(obj);
    v = flipud(v);
    if video_type == 1
        image_offset = size(v,2);
    elseif video_type == 2
        image_offset = 0;
    end
    
    % load correct Inscopix movie frame
    try
        frame = double(h5read(infile,'/Object',[1 1 ceil(t(i)*aviobj.FrameRate)+FToffset 1],[Xdim Ydim 1 1]));
    catch
        return
    end
    % rescale the frame to [0 256]
    frame(find(frame < clims(1))) = clims(1);
    frame(find(frame > clims(2))) = clims(2);
    frame = frame-clims(1);
    %keyboard;
    frame = uint8(frame/(clims(2)-clims(1))*256);
    
    if video_type == 1
        combXdim = max(size(v,1),size(frame,1));
        combYdim = image_offset+size(frame,2);
    elseif video_type == 2
        combXdim = size(frame,1);
        combYdim = size(frame,2);
    elseif video_type == 3
        combXdim = size(v,1);
        combYdim = size(v,2);
    end
    
    temp = zeros(combXdim,combYdim,3);
    if video_type ~= 2 % Plot tracking
        temp(1:size(v,1),1:size(v,2),:) = v;
    end
    
    if video_type ~=3 % Plot brain imaging
        temp(1:size(frame,1),image_offset+1:end,3) = frame;
        temp(1:size(frame,1),image_offset+1:end,2) = frame;
        temp(1:size(frame,1),image_offset+1:end,1) = frame;
    end
    image(uint8(temp));
    axis equal;axis tight; axis off;hold on;
    Xa = get(gca,'XLim');
    Ya = get(gca,'YLim');
    
    if video_type == 1
        set(gcf,'Position',[1          41        1920         964]);
    else
        set(gcf,'Position',[1          41        1050        964]);
    end
    % plot trajectory, hold on
    %     plot(xAVI,yAVI,'-','Color',[0.2 0.2 0.2]);hold on;axis tight;
    
    % find active neurons
    an = find(FT(:,i));
    hold on;
    if video_type == 1 || video_type == 2 % Plot Neuron outlines
        for j = 1:length(an)
            plot(xOutline{an(j)}+image_offset,yOutline{an(j)},'-r','LineWidth',3,'Color',colors(an(j),:));
        end
    end
    
    xg = [];
    yg = [];
    
    if video_type ~=2
        % for each active neuron
        for j = an'
            
            if ((pval(j) > 0.95) & (nt(j) >= 3))
                
                plot(xt{j},yt{j},'Color',colors(j,:),'LineWidth',5);
                
            end
        end
        end
    
    if video_type ~= 2
        % plot mouse marker
        if (isrunning(i))
            mf = 'r';
        else
            mf = 'k';
        end
        
        plot(Xb2AVI(Xbin(i)),Yb2AVI(Ybin(i)),'ok','MarkerSize',10,'MarkerFaceColor',mf);
        set(gca,'XLim',Xa,'YLim',Ya);
    end
    
    
    % getframe
    F = getframe(gcf);
    % write to avi
    writeVideo(aviobj,F);
    hold off;
    gcf;
    
end

end
