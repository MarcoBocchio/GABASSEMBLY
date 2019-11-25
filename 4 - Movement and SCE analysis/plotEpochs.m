function [] = plotEpochs(logicalArray,time);
hold all;

if class(logicalArray) == 'double'
    logicalArray = logical(logicalArray);
end

if nargin < 2
    time = 1:length(logicalArray);
end

    for trialCount = 1:size(logicalArray,1)
        epochPos = time(logicalArray(trialCount, :));
        for epochCount = 1:length(epochPos)
            plot([epochPos(epochCount) epochPos(epochCount)], ...
                [trialCount-0.4 trialCount+0.4], '-','k'); 
        end
    end
    
xlim([0 max(time)]);

if nargin < 2
    xlabel('Time (ms)');
else
    xlabel('Frame #')
end

