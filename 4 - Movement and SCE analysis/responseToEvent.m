function [spikeHist_edges,spikeProb,spikeZscore,signPosBins,signNegBins,excitResponse,inhibResponse] = responseToEvent (spikeCount,eventTimes,PSTHwindow,nBins,si, option, plotting, offset)

%% >>> OPERATION >>>
% Define responsive and non responsive cells with z-score analysis of PSTH.

% >>> INPUTS >>>
% NAME             TYPE, DEFAULT        DESCRIPTION
% spikeCount       double               spike raster array
% eventTimes       double               reference event onset (stimulation,behavioural epoch)
% PSTHwindow       scalar, 50           size of window around reference event (in data points)
% nBins            scalar, 30           number of bins for histogram
% si               scalar,1             sampling interval (in ms). Set as 1 if not provided
% option           case
%                   'post'              finds response after stimulus
%                   'during'            finds response around stimulus
% offset           scalar, 0            offset for baseline (e.g. set on -2
%                                       to search response starting from
%                                       second bin before the stimulus)


% >>> OUTPUTS >>>
% NAME             TYPE                 DESCRIPTION
% spikeProb        double               spike probability histogram
% spikeZscore      double               spike z-score histogram
% signBins         double               location of significant bins
% response         logical              outcome (1: response, 2: no response)

% 
% Marco Bocchio, 16/7/1

%%
% set default values
if nargin < 8
    offset = 0;
    if nargin < 7
        plotting = 'plotting';
        if nargin < 6
            option = 'post';  
            if nargin < 5
                si = 1;
                    if nargin < 2
                    PSTHwindow = 50;
                    nBins = 30;
                    end
            end
        end
    end
end



% number of trials (i.e. events or stimulations)
trialN = length(eventTimes);

% window around event or stimulation onset
%lags = linspace(-PSTHwindow,PSTHwindow,(PSTHwindow*2)+1);

%pre-allocation
alignedSpikeCount = zeros(trialN,PSTHwindow*2+1);

% spike count sorted by trial (each trial is an event or stimulation onset)
for i = 1:trialN
    eventCounter = eventTimes(i);
    if eventCounter <= PSTHwindow
        continue
    end
    if (eventCounter+PSTHwindow) >= size(spikeCount,2)
        continue
    end
    trial = spikeCount(eventCounter-PSTHwindow:eventCounter+PSTHwindow);
    alignedSpikeCount(i,:) = trial;
end

% spike times aligned at event onset
[alignedSpikeTimes] = rasterToTimes(alignedSpikeCount,si,PSTHwindow+1);

% concatenated spike times
alignedSpikeTimes=cell2mat(alignedSpikeTimes');

% PSTH of spike number / probability
spikeHist_edges = linspace(-PSTHwindow*si,PSTHwindow*si,nBins);
[spikeHist] = histcounts(alignedSpikeTimes,spikeHist_edges);
spikeProb = spikeHist./trialN;
baselineBins = find(spikeHist_edges(2:end) < 0); %define baseline length
baselineBins = baselineBins(1:end+offset);

if mean(spikeProb(baselineBins)) == 0 %set spike probability in baseline to a very small number (it can't be exactly zero for z-score calculation)
    spikeProb(baselineBins) = 0.001; % default: 0.001
end

spikeZscore = (spikeProb-(mean(spikeProb(baselineBins(1:end-2)))))./(std(spikeProb(baselineBins(1:end-2)))); % z-scored PSTH


switch option
    case 'post' % a significant response is searched after the stimulus
        [signPosBins] = consecAboveThresh(spikeZscore(max(baselineBins):end),2,2); %significant consecutive bins are searched after stimulus onset (excitatory response)
        [signNegBins] = consecAboveThresh(-spikeZscore(max(baselineBins):end),2,2); %significant consecutive bins are searched after stimulus onset (inhibitory response)
        consecBinsThresh = 1; %threshold of consecutive bins to define response
        
    case 'during' %a significant response locked to the stimulus (i.e. SCEs) is searched 
        
        if max(spikeHist_edges)<(2*si) % avoids error in case firing rate of the cell is too low (automatically sets the cell as non responsive)
        excitResponse = logical(false);
        inhibResponse = logical(false);
        signPosBins = [];
        signNegBins = [];
        return
        end

        [signPosBins] = find(spikeZscore(max(baselineBins):max(baselineBins)+1)>3); % one significant bin is searched in the 3 bins around the stimulus onset (excitatory response)
        [signNegBins] = find(-spikeZscore(max(baselineBins):max(baselineBins)+1)>3); % one significant bin is searched in the 3 bins around the stimulus onset (inhibitory response)
        consecBinsThresh = 0; %threshold of consecutive bins to define response
end

% define excitation
if length(signPosBins)>consecBinsThresh
    excitResponse = logical(true);
else
    excitResponse = logical(false);
end


% define inhibition
if length(signNegBins)>consecBinsThresh
    inhibResponse = logical(true);
else
    inhibResponse = logical(false);
end

switch plotting
    case 'plotting'
        
        subplot(2,1,1)
        bar(spikeHist_edges(2:end),spikeProb);
        ylabel('Spike probability');
        xlabel('Time (ms)')
        subplot(2,1,2)
        bar(spikeHist_edges(2:end),spikeZscore);
        ylabel('Z-score');
    case 'noPlotting'
    otherwise  
end






end