function [ output_args ] = MPFM_mousevid(out_avifile,Pix2Cm)
% make a fun movie of our mouse and his placefields
%   Detailed explanation goes here

close all;

if (~exist('Pix2Cm'))
    Pix2Cm = 0.15;
    display('assuming room 201b');
    % factor for 201a is 0.0709
else
    if (strcmp(RoomStr,'201a'))
        Pix2Cm = 0.0709;
        display('Room 201a');
    end
end

aviSR = 30.0003;

try
    h1 = implay('Raw.AVI');
    obj = VideoReader('Raw.AVI');
catch
    avi_filepath = ls('*.avi');
    h1 = implay(avi_filepath);
    disp(['Using ' avi_filepath ])
    obj = VideoReader(avi_filepath);
end

load PlaceMaps.mat;
load PFstats.mat;

NumFrames = length(x);
NumNeurons = length(NeuronImage);

figure;
set(gcf,'Position',[534 72 1171 921]);

aviobj = VideoWriter(out_avifile);
aviobj.FrameRate = 20;
aviobj.Quality = 10;
open(aviobj);

% assign each neuron a color
colors = rand(NumNeurons,3);

% convert x and y to .avi coords
xAVI = x/Pix2Cm*0.625;
yAVI = y/Pix2Cm*0.625;

% convert Xbin and Ybin to x and y
Xd = Xedges(2)-Xedges(1);
Yd = Yedges(2)-Yedges(1);

for i = 1:length(Xedges)
    Xb2AVI(i) = (Xedges(i)+Xd/2)/Pix2Cm*0.625;
end

for i = 1:length(Yedges)
    Yb2AVI(i) = (Yedges(i)+Yd/2)/Pix2Cm*0.625;
end

for i = 1:NumFrames
    
    % load correct movie frame
    % calculate correct frame based on iteration and offsets
    obj.currentTime = (i/20+(FToffset-16)/20);
    v = readFrame(obj);
    
    imagesc(flipud(v));axis equal;hold on;
    
    
    
    set(gcf,'Position',[534 72 1171 921]);
    % plot trajectory, hold on
%     plot(xAVI,yAVI,'-','Color',[0.2 0.2 0.2]);hold on;axis tight;
     Xa = get(gca,'XLim');
     Ya = get(gca,'YLim');
    
    
    
    % find active neurons
    an = find(FT(:,i));
    
    % for each active neuron
    for j = an'
        % get PF outline (if avail)
        WhichField = MaxPF(j);
        temp = zeros(size(TMap{1}));
        tp = PFpixels{j,WhichField};
        try
        temp(tp) = 1;
        catch
            keyboard;
        end
        % plot PF outline (using correct color)
        b = bwboundaries(temp,4);
        if (~isempty(b))
            
            yt = Yb2AVI(b{1}(:,2));
            xt = Xb2AVI(b{1}(:,1));
            xt= xt+(rand(size(xt))-0.5)/2;
            yt= yt+(rand(size(yt))-0.5)/2;
            %colors(j,:)
            plot(xt,yt,'Color',colors(j,:),'LineWidth',5);
            
            
        end
    end
    
    % plot mouse marker
    if (isrunning(i))
        mf = 'r';
    else
        mf = 'k';
    end
    
    plot(Xb2AVI(Xbin(i)),Yb2AVI(Ybin(i)),'ok','MarkerSize',10,'MarkerFaceColor',mf);
    set(gca,'XLim',Xa,'YLim',Ya);
    
    % getframe
    F = getframe(gcf);
    % write to avi
    writeVideo(aviobj,F);
    hold off;
    gcf;
end
