function [spikeTimes,spikeTimes_conc] = rasterToTimes(spikenums,si,offset)
%% >>> OPERATION >>>
% Extract spike times from a logical spike raster (in which each row is a cell) 

% >>> INPUTS >>>
% NAME             TYPE, DEFAULT        DESCRIPTION
% spikenums        double               spike raster matrix 
% si               scalar,1             sampling interval (in ms). Set as 1 if not provided
% offset           scalar,0             offset to zero spike times in
%                                       correspondence of an event. Set as zero if not provided

% >>> OUTPUTS >>>
% NAME             TYPE                 DESCRIPTION
% spikeTimes       cell                 spike times of each cell
% spikeTimes_conc  double               concatenate spike times
%                

% 
% Marco Bocchio, 10/7/19

if nargin < 3

    if exist('si') == 0
        si = 1;
    end

    if exist('offset') == 0
        offset = 0;
    end
end

rowN = size(spikenums,1);

%x/y coordinates of spikes in 2D array
[spikeTimes_rowN,spikeTimes_times]=find(spikenums==1);

% convert to actual spike 'times' if sampling interval is provided
if nargin > 1
    spikeTimes_times = spikeTimes_times * si;
end

%store spike times of individual cells in cell array
for rowCounter = 1:rowN
    spikeTimes_rowIndex = spikeTimes_times(spikeTimes_rowN == rowCounter);
    spikeTimes_rowIndex = spikeTimes_rowIndex - (offset*si);
    spikeTimes{rowCounter,:} = spikeTimes_rowIndex';
end

end

