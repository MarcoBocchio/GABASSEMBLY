function [spikenums_shuffled,syncEvents_shuffled,SCEthresh] = circshiftRaster (spikenums,number_of_surrogates,plotting);

%% >>> OPERATION >>>
% Reshuffles spike raster to generate:
% - surrogate spike raster (same size as original raster)
% - surrogates for synchronous events (n surrogates x n data points 2D array)

% If surrogates for synchronous events are needed, input spike raster needs
% be a logical array (and not a spike count array)

%% >>> INPUTS >>>
% NAME                  TYPE, DEFAULT        DESCRIPTION
% spikenums             double               spike raster matrix 
% number_of_surrogates  scalar, 1000         
% 

%% >>> OUTPUTS >>>
% NAME                  TYPE, DEFAULT        DESCRIPTION
% spikenums_shuffled    double               shuffled spike raster 2D array 
% syncEvents_shuffled   double               shuffled synchronous events 2D array 
% % 
%% 
% Marco Bocchio, 8/7/19

if nargin < 3
    plotting = 'noPlotting';
end

if nargin < 2
    number_of_surrogates = 1000;
end
         
    for surr_counter = 1:number_of_surrogates
        spikenums_shuffled=zeros(size(spikenums)); %pre-allocation
        for neuron_counter = 1:size(spikenums,1)
            spikenums_shuffled(neuron_counter,:)=circshift(spikenums(neuron_counter,:), randi([1, size(spikenums,2)]));
        end
        syncEvents_shuffled_idx=sum(spikenums_shuffled,1); % synchronous events for each cell in the loop
        syncEvents_shuffled{surr_counter}=syncEvents_shuffled_idx; % synchronous events (n surrogates x n data points 2D array)
    end
    syncEvents_shuffled = cell2mat(syncEvents_shuffled); %concatenation of 1000 reshuffled syncEvents (sum of activity along columns) arrays
    SCEthresh=prctile(syncEvents_shuffled,99); %99th percentile of reshuffled distribution used as threshold
        
switch plotting
    case 'plotting'
        figure
        subplot(2,1,1)
        plotSpikeRaster(logical(spikenums),'PlotType','vertline');
        title('Original raster');
        subplot(2,1,2)
        plotSpikeRaster(logical(spikenums_shuffled),'PlotType','vertline');
        title('Reshuffled raster');
    case 'noPlotting'
    otherwise
end


    

end
    
