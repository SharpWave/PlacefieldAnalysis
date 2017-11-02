function [DaveTMap,DaveTMapRaw] = DavePF(FT,PositionVector,isrunning,NumXbin,NumYbin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NumPos = NumXbin * NumYbin;

NumNeurons = size(FT,1);
NumSamples = size(FT,2);
LastPV = 0;
Curr = 0;

% first, make NewFT and new PV

for t = 1:NumSamples
    if(~isrunning(t))
      continue;
    end

    if ((LastPV == 0) || (PositionVector(t) ~= LastPV))
        % new position bin
        Curr = Curr + 1;
        CurrT(Curr) = t;
        
        LastPV = PositionVector(t);
        NewPV(Curr) = LastPV;
        NewFT(1:NumNeurons,Curr) = 0;
    end
    
    for i = 1:NumNeurons
        if(NewFT(i,Curr) == 0)
            NewFT(i,Curr) = FT(i,t);
        end
    end
end

NumVisits = Curr;
TrCount = zeros(NumNeurons,NumPos);
VisitCount = zeros(1,NumPos);

for t = 1:NumVisits
    VisitCount(NewPV(t)) = VisitCount(NewPV(t)) + 1;
    for i = 1:NumNeurons
        TrCount(i,NewPV(t)) = TrCount(i,NewPV(t)) + NewFT(i,t);
    end
end


sm_gauss = fspecial('gaussian',[round(8*5,0), round(8*5,0)],5);

NewPosMap = reshape(VisitCount,[NumXbin,NumYbin]);
NewPosMap = ordfilt2(NewPosMap,100,ones(10));
NewPosMap = imfilter(NewPosMap,sm_gauss);

for i = 1:NumNeurons
    DaveTMapRaw{i} = reshape(TrCount(i,:),[266,169]);
    DaveTMap{i} = ordfilt2(DaveTMapRaw{i},100,ones(10,10));
    DaveTMap{i} = imfilter(DaveTMap{i},sm_gauss);
    DaveTMap{i} = DaveTMap{i}./NewPosMap;
    
end    


