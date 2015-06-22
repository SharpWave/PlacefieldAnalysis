function [ ] = batch_align_pos(base_struct, reg_struct)
% Aligns position data so that every session has the same bounds on the
% occupancy map and we can easily do correlations and other comparisons
% between firing across sessions
%
% INPUTS: mirror MD from MakeMouseSessionList, but must include at least
% .Animal, .Date, .Session, AND .Room fields

%% 1: Load all sessions, and align to imaging data

% Dumpt everything into one structure for future ease
sesh(1) = base_struct;
sesh(2:length(reg_struct) + 1) = reg_struct;

currdir = cd;
for j = 1: length(sesh)
    ChangeDirectory(sesh(j).Animal, sesh(j).Date ,sesh(j).Session);
    if ~isempty(regexpi(sesh(1).Room,'201b'))
        Pix2Cm = 0.15;
        disp(['Using 0.15 for Pix2Cm for ' sesh(j).Date ' Session ' num2str(sesh.Session)])
    else
        Pix2Cm = [];
        disp('Need room to get Pix2Cm')
    end
    load('ProcOut.mat', 'FT')
    [x,y,~,FTalign,FToffset,FToffsetRear] = AlignImagingToTracking(Pix2Cm,FT);
    sesh(j).x = x;
    sesh(j).y = y;
end

keyboard
%% 2: Align position data for each session to the base session by using the 95% occupancy limits, save as Pos_align.mat
% Include base session in Pos_align for future reference

%% 4: Concatenate ALL position data into one X and one Y vector, and get Xedges and Yedges based on this

%% 5: Save Xedges in base session for future reference along with all sessions aligned to it.


end

