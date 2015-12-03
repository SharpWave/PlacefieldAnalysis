function [] = PostHocCompare(h5file,CellsToBrowse,varargin)
% compare cell outlines obtained via regular Tenaspis vs.
% PostTenaspisCorrections.m

load PlaceMaps.mat;
load ProcOut.mat;
load Post_FT.mat;

c = colormap;
c = [[0 0 0];c];

NumNeurons = length(NeuronImage);
NumFrames = size(FT,2);

if (nargin < 2)
    CellsToBrowse = 1:NumNeurons;
end

for i = CellsToBrowse
    i
    % Plot # 1 : FT Movie frames
    ax1 = subplot(1,3,1);
    activeframes = find(FT(i,:) == 1);
    avgframe{i} = zeros(size(NeuronImage{1}));
    for j = activeframes

            
        avgframe{i} = avgframe{i} + double(loadframe(h5file,j));
    end
    
    ae = NP_FindSupraThresholdEpochs(FT(i,:),eps);
    NumTransients = size(ae,1);
    avgframe{i} = avgframe{i}./length(activeframes);
    
    imagesc(avgframe{i});hold on;caxis([-500 500]);colormap(ax1,gray);plot(xOutline{i},yOutline{i},'-r');
    title([int2str(NumTransients),' transients']);hold off;
    

    
    % plot 2 - PFT Movie frames
    ax2 = subplot(1,3,2);
    activePframes = find(PFT(i,:) == 1);
    avgPframe{i} = zeros(size(NeuronImage{1}));
    for j = activePframes
         
        avgPframe{i} = avgPframe{i} + double(loadframe(h5file,j));
    end
    
    ae = NP_FindSupraThresholdEpochs(PFT(i,:),eps);
    NumTransients = size(ae,1);
    avgPframe{i} = avgPframe{i}./length(activePframes);
    
    imagesc(avgPframe{i});hold on;caxis([-500 500]);colormap(ax2,gray);plot(xOutline{i},yOutline{i},'-r');
    title([int2str(NumTransients),' transients']);hold off;
    
    ax3 = subplot(1,3,3);
    imagesc(avgPframe{i}-avgframe{i});caxis([-400 400]);colormap(ax3,'default');colorbar;hold on;plot(xOutline{i},yOutline{i},'-r');hold off;
    pause;
end

