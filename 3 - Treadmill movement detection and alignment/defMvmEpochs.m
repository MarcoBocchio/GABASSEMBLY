function [mvmEpochsLogic, mvmEpochsIndex, mvmOnsetIndex, mvmOffsetIndex, restEpochsLogic, restEpochsIndex, time] = defMvmEpochs (mvmTrimmed,spikenums,si_abf,si_img,spanThreshold_1,spanThreshold_2,rasterOptions)

% >>> OPERATION >>>
% Use wheel sensor data to define epochs with movement
% Only periods with flat signal > spanThreshold2 (usually 200 ms) are selected as
% 'rest' epochs. Shorter flat periods intermingled with movement periods
% are included in nearby movement periods
%
% % >>> INPUT VARIABLES >>>
% NAME             TYPE, DEFAULT        DESCRIPTION
% spikenums        double               spike raster matrix 
% si_img           scalar, 50.75        sampling interval for imaging data
% si_abf           scalar, 0.2          sampling interval of abf file         
% spanThreshold_1  scalar, 1000         initial threshold to find
%                                       consecutive periods of rest
% spanThreshold_2  scalar, 200          second threshold to label rest
%                                       periods as such only if their
%                                       duration is > spanThreshold_2
%%
% Marco Bocchio, 4/7/19

spanThreshold_1 = spanThreshold_1 / si_abf;
spanThreshold_2 = spanThreshold_2 / si_abf;


%% set consecutive ones between movement bouts as zero (set all rest periods as zero)
[restIndex_1] = consecAboveThresh(mvmTrimmed,0.05,spanThreshold_1);
mvmTrimmed(restIndex_1) = 0; 
mvmTrimmed = mvmTrimmed > 0.05;


%% set short epochs of zeros as ones (include in nearby movement epochs)
rest=abs(mvmTrimmed-1); % vector in which ones define rest periods
[restIndex_2] = consecAboveThresh(rest,0.05,spanThreshold_2); %find consecutive rest periods (ones) above a threshold
mvmTrimmed=ones(length(mvmTrimmed),1);
mvmTrimmed(restIndex_2)=0;

%% downsample movement vector to match spike raster data points
mvmEpochsLogic = round(resample(mvmTrimmed,length(spikenums),length(mvmTrimmed)));
if length(mvmEpochsLogic)>length(spikenums)
    mvmEpochsLogic = mvmEpochsLogic(1:end-(length(mvmEpochsLogic)-length(spikenums)));
end
mvmEpochsIndex = find (mvmEpochsLogic == 1);

%% movement onset and offset
mvmOnsetIndex = (find(diff(mvmEpochsLogic)==1)+1);
mvmOffsetIndex = (find(diff(mvmEpochsLogic)==-1)+1);

%% rest epochs
restEpochsLogic = ~mvmEpochsLogic;
restEpochsIndex = find (restEpochsLogic == 1);

%% calculate time vector
time = 0:si_img*10^-3:(length(spikenums)-1)*si_img*10^-3; % in s

%% plotting
if rasterOptions.epochs.epochs2 == 'mvm';
    rasterOptions.epochs.epochs2 = mvmEpochsIndex;
end


figure
plotRaster(spikenums,rasterOptions);
ylabel('Cell #');
title('Spike raster and movement epochs');


end