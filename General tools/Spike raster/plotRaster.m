function [] = plotRaster(spikenums,rasterOptions);

%% >>> OPERATION >>>
%   Plots spike raster for a 2D array in which rows are cells and columns
%   are data points.
%   Options:    - plot time vector as x
%               - highlight specific cells in red
%               - highlight specific time windows (e.g. behavioural epochs,
%                   stimuli) with coloured vertical lines behind the raster

% >>> INPUTS >>>
% NAME             TYPE, DEFAULT        DESCRIPTION
% spikenums        double               spike raster matrix 
% options          struct               set parameters (see below
% 

% >>> OPTIONS >>>
% NAME                      DESCRIPTION
% rasterOptions.time                 time vector
% rasterOptions.cellColour           number of cell(s) to be highlighted in red
% rasterOptions.epochs               period(s) of time to be highlighted (3 max)
% 
% rasterOptions.epochs to be provided as:
% rasterOptions.epochs.epochs1
% rasterOptions.epochs.epochs2
% rasterOptions.epochs.epochs3   

%  
%
% Marco Bocchio, 10/7/19

%%
hold all;

% default options if rasterOptions not provided
if nargin < 2
    rasterOptions.time = [];    
    rasterOptions.cellColour = [];
    rasterOptions.epochs.epochs1 = [];
end

% default options if not all options are provided inside rasterOptions
if isfield(rasterOptions,'time') == 0; 
    rasterOptions.time = [];  
end

if isfield(rasterOptions,'cellColour') == 0; 
    rasterOptions.cellColour = [];  
end

if isfield(rasterOptions,'epochs') == 0; 
    rasterOptions.epochs.epochs1 = [];
    rasterOptions.epochs.epochs2 = [];
    rasterOptions.epochs.epochs3 = [];
end


% make spike matrix logical if not already
if class(spikenums) == 'double'
    spikenums = logical(spikenums);
end

% create  time vector if not existent (data points will be plotted instead
% of time points)
if isempty(rasterOptions.time) == 1;
    rasterOptions.time = 1:length(spikenums);
end

% highlight desired cell(s) in red
if isempty(rasterOptions.cellColour) == 0;
    for trialCounter = 1:size(spikenums,1)
        spikePos = rasterOptions.time(spikenums(trialCounter, :));
        trialToColourCode = any(rasterOptions.cellColour(:) == trialCounter); %check if current trial in the loop is among the ones to colour code
        for spikeCount = 1:length(spikePos)
            if trialToColourCode == true %if current trial is to colour code
                plot([spikePos(spikeCount) spikePos(spikeCount)], ...
                [trialCounter-0.4 trialCounter+0.4], 'r'); %colour in red
            else
                plot([spikePos(spikeCount) spikePos(spikeCount)], ...
                [trialCounter-0.4 trialCounter+0.4], 'k'); %else in black
            end
        end
    end
    
else
    
   for trialCounter = 1:size(spikenums,1)
        spikePos = rasterOptions.time(spikenums(trialCounter, :));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCounter-0.4 trialCounter+0.4], 'k'); %else in black
        end
   end 
    
end
    
ylim([0 size(spikenums, 1)+1]);
xlim([0 max(rasterOptions.time)]);

% label x axis depending on whether a time vector was provided as input or
% not
if max(rasterOptions.time) == length(spikenums);
    xlabel('Frame #')
else
    xlabel('Time (ms)');
end

% highlight specific epochs
%if isempty(rasterOptions.epochs.epochs1) == 0;    
    if isfield(rasterOptions.epochs,'epochs1') == 1;
       for i = 1:length(rasterOptions.epochs.epochs1)
          line([rasterOptions.epochs.epochs1(i) rasterOptions.epochs.epochs1(i)],[0.5 size(spikenums,1)+0.5],'Color','green');
       end
    end
    
     if isfield(rasterOptions.epochs,'epochs2') == 1;
       for j = 1:length(rasterOptions.epochs.epochs2)
          line([rasterOptions.epochs.epochs2(j) rasterOptions.epochs.epochs2(j)],[0.5 size(spikenums,1)+0.5],'Color',[0.8 0.8 0.8]);
       end
     end
    
    if isfield(rasterOptions.epochs,'epochs3') == 1;
       for k = 1:length(rasterOptions.epochs.epochs3)
          line([rasterOptions.epochs.epochs3(k) rasterOptions.epochs.epochs3(k)],[0.5 size(spikenums,1)+0.5],'Color','yellow');
       end
    end
        
    
    set(gca,'children',flipud(get(gca,'children')))
    ylabel('Cell #');
    
%end
    
end