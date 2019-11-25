function [alignedSpikeCount] = sortByTrial (spikeCount,eventTimes,PSTHwindow,offset)

%% >>> OPERATION >>>
% Sort 2D spike matrix into 3D spike matrix based on event index (e.g.
% locomotion onset). Trials are arranged in 3rd dimension.


% Marco Bocchio, 22/11/19

%%

% number of trials (i.e. events or stimulations)
trialN = length(eventTimes);

% window around event or stimulation onset
%lags = linspace(-PSTHwindow,PSTHwindow,(PSTHwindow*2)+1);

%pre-allocation
alignedSpikeCount = zeros(size(spikeCount,1),PSTHwindow-offset+1,trialN);

% spike count sorted by trial (each trial is an event or stimulation onset)
for trialCounter = 1:trialN
    eventCounter = eventTimes(trialCounter);
    if eventCounter <= PSTHwindow
        continue
    end
    if (eventCounter+PSTHwindow) >= size(spikeCount,2)
        continue
    end
    alignedSpikeCount(:,:,trialCounter) = spikeCount(:,eventCounter+offset:eventCounter+PSTHwindow);
end

end 

