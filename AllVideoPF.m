function [ ] = VideoPF(NeuronIdx)

close all;
load PlaceMaps.mat;

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

NumFrames = length(x);
NumNeurons = length(NeuronImage);
Xdim = size(NeuronImage{1},1);
Ydim = size(NeuronImage{1},2);

tempv = zeros(obj.Height,obj.Width,3,'double');
allv = zeros(obj.Height,obj.Width,3,'double');
NumUsed = 0;
for i = 1:NumFrames
    
    % load correct Plexon movie frame
    % calculate correct frame based on iteration and offsets
    obj.currentTime = aviFrame(i);
    v = readFrame(obj);
    v = flipud(v);
    allv = allv+double(v);
 
end

allv = allv./NumFrames;

save allv.mat allv;

end
