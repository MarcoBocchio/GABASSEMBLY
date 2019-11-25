function [binnedRaster] = binRaster(rasterToBin, binEdges);

%%% >>> OPERATION >>>
% Binning of spike raster (moving sum of n nearby points, where
% n=binEdges).
% Probably slower than using sepblockfun function, but works even when
% binEdges is not multiple of the number of time points.

% >>> INPUT VARIABLES >>>
% NAME             TYPE, DEFAULT        DESCRIPTION
% rasterToBin      double               spike raster matrix (cells along rows and times/frames along columns)
% binEdges         scalar,              binning (i.e. nearby points to be summed)
% 

%% 
% Marco Bocchio, 10/7/19

cellN = size(rasterToBin,1);
frameN = size(rasterToBin,2);
newFrameN = floor(frameN  / binEdges);

binnedRaster = zeros(cellN,newFrameN);

for i = 1:1:cellN
    binaryArrayToBin = rasterToBin(i,:);
    binnedArray_i = sum(reshape(binaryArrayToBin(1:binEdges * floor(numel(binaryArrayToBin) / binEdges)), binEdges, []),1);
    binnedRaster(i,1:length(binnedArray_i)) = binnedArray_i;
end

end