function [ output_args ] = ComparePFBrowse(animal_id,date1,sess1,date2,sess2,minpval)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ChangeDirectory(animal_id,date2,sess2);
load PlaceMaps.mat;

CellsToTry = find(pval >= minpval);


end

