function [nCellsSCE,SCEindex,SCE_logic,nSpikes_SCE,SCEthresh,spikenums_binned] = defSCEs(spikenums,si_img,minDistSCEs,binning,rasterOptions)
%% >>> OPERATION >>>
% Defines SCEs in spike raster. 1000 surrogates matrix of synchronous events
% is generated. The 95th/99th percentile of each column (data point) is
% defined as a dynamic threshold to detect SCEs.

%% >>> INPUT VARIABLES >>>
% NAME             TYPE, DEFAULT        DESCRIPTION
% spikenums        double               spike raster matrix 
% si_img           scalar, 101.5        sampling interval for imaging recording (frame rate, in ms)
% minDistSCEs      scalar, 1000         minimum distance between SCEs (in ms)
% binning          scalar, 2            bin edges to "downsample" spike raster before SCE detection

%% 
% Marco Bocchio, 10/7/19

si_SCE = si_img * binning;



% Bin raster
[spikenums_binned] = binRaster(spikenums, binning); %binned spike raster
spikenums_binned_logic = spikenums_binned;
spikenums_binned_logic(spikenums_binned_logic > 1) = 1; %remove multiple spikes per frame to calculate co-active cells

time = 0:1:size(spikenums_binned,2)-1;
time = time * (si_SCE*10^-3);


% Surrogate raster
[spikenums_shuffled,syncEvents_shuffled,SCEthresh] = circshiftRaster (spikenums_binned_logic,1000);

% Synchronous events
syncEvents = sum(spikenums_binned_logic,1);

%SCEs in spike raster (plot only)
figure;
subplot(212)
findpeaks(syncEvents,'MinPeakDistance',round(minDistSCEs / si_SCE),'MinPeakHeight',SCEthresh);
title('Detected SCEs');

% SCEs in spike raster (save peaks and indexes)
[nCellsSCE,SCEindex] = findpeaks(syncEvents,'MinPeakDistance',round(minDistSCEs / si_SCE),'MinPeakHeight',SCEthresh);
SCE_logic = zeros(1,size(spikenums_binned_logic,2));
SCE_logic(SCEindex) = 1;

% number of spikes / SCE
nSpikes_SCE = sum(spikenums_binned(:,SCEindex),2)./length(SCEindex);

if rasterOptions.epochs.epochs1 == 'SCE';
    rasterOptions.epochs.epochs1 = SCEindex;
end

SCEfreq = length(SCEindex)/max(time);
disp(['SCE frequency:' num2str(round(SCEfreq,2)) ' Hz']);

% Plotting
subplot(211)
plotRaster(spikenums_binned_logic,rasterOptions);
ylim([0 size(spikenums,1)])
ylabel('Cell #');




end

