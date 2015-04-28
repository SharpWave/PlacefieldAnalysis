function [] = PFA(roomstr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Step 1: Calculate Placefields
CalculatePlacefields(roomstr);

% Step 2: Calculate Placefield stats
PFstats();

% Step 3: Extract Traces
ExtractTracesProc('D1Movie.h5',infile)

% Step 4: Browse placefields
PFbrowse('D1Movie.h5');

end

