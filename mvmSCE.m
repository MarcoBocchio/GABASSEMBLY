%--------------------------------------------------------------------------
% Pipeline to analyse single cell activities during movement and SCEs
%--------------------------------------------------------------------------
%
% >>> STRUCTURE >>>
%   -  Load ABF file to extract treadmill movement and 'imaging on' periods.
%       Exclude treadmill data during 'imaging off' periods
%   -  Define movement and rest epochs
%   -  Measure firing rates during movement and rest
%   -  Define movement modulated cells based on z-score and spike ratio
%       analyses
%   -  Detect Synchronous Calcium Events (SCEs) and SCE firing modulation

% >>> REQUIRED THIRD-PARTY CODES >>>
%   -   abfload (F. Collman). https://github.com/fcollman/abfload

% >>> TO DO LIST >>>
%   -  Treadmill position and speed analysis
%   -  Spatial selectivity (i.e. place cells) analysis
%   -  Analysis of LFP data


%--------------------------------------------------------------------------
% Marco Bocchio, updated 22/11/2019
%--------------------------------------------------------------------------


%% Load and trim ABF

imagingCh = 1; %imaging channel in abf file
mvmCh = 4;     %locomotion channel in abf file

if isstruct(spikenums)==0; %if spikenums is not a struct array yet, convert it to struct
    spikenums2 = spikenums;
    clear spikenums;
    spikenums.tot.orig = spikenums2;
    clear spikenums2;
end


[imagingTrimmed, mvmTrimmed, fs_abf, si_abf] = loadMvm(spikenums.tot.orig,imagingCh,mvmCh);

if exist('fileName') == 0;
    fileName = input('Type filename:','s');
end

if exist('si_img') == 0;
    si_img = input('Type sampling interval for imaging:','s');
end

si_img = si_img/2; %double resolution after MCMC

fileName_mvm=strcat(fileName,'_mvm.mat');
save(fileName_mvm,'fileName', 'mvmTrimmed','fs_abf','si_abf','si_img','-v7.3');
%clear fileName_mvm imagingTrimmed timeStart timeEnd;

%% Define movement and rest epochs and align them to spike raster
clear rasterOptions; % options for plotting the raster
rasterOptions.cellColour = nCells.class2; % cells to highlight in binned spike raster
rasterOptions.epochs.epochs2 = 'mvm'; % highlight movement in raster

[mvmEpochsLogic, mvmEpochsIndex, mvmOnsetIndex, mvmOffsetIndex, restEpochsLogic, restEpochsIndex, time] = defMvmEpochs (-mvmTrimmed,spikenums.tot.orig,si_abf,si_img,1000,1000,rasterOptions);

if exist('fileName') == 0;
    fileName = input('Type filename:','s');
end

% remove pre-movement period (20 frames) from rest periods (avoid
% contamination for SCE and assembly analysis)
for mvmIndexCounter = 1:length(mvmOnsetIndex)
preMvmIndexTemp = mvmOnsetIndex(mvmIndexCounter)-20:mvmOnsetIndex(mvmIndexCounter)-1;
    if mvmIndexCounter == 1
        preMvmIndex = preMvmIndexTemp;
    else
        preMvmIndex = [preMvmIndex preMvmIndexTemp];
    end
end

clear preMvmIndexTemp

restEpochsOrigLogic = restEpochsLogic; % save a copy of original rest time points
restEpochsOrigIndex = restEpochsIndex;

restEpochsLogic(preMvmIndex)=[]; % remove pre-movement periods from rest
restEpochsIndex = find (restEpochsLogic == 1);

fileName_mvmEpochs=strcat(fileName,'_mvmEpochs.mat');
save(fileName_mvmEpochs,'fileName', 'mvmEpochsLogic','mvmEpochsIndex','mvmOnsetIndex','mvmOffsetIndex', 'restEpochsLogic', 'restEpochsIndex','si_img','time','-v7.3');
%clear fileName_mvmEpochs mvmTrimmed fs_abf si_abf;

%% Rates during rest and movement
% spike rates during rest and movement
nSpikes.tot = sum(spikenums.tot.orig,2); % n of spikes
nSpikes.rest = sum(spikenums.tot.orig(:,restEpochsIndex),2);
nSpikes.mvm = sum(spikenums.tot.orig(:,mvmEpochsIndex),2);

spikeRate.tot = nSpikes.tot./(size(spikenums.tot.orig,2)*si_img*10^-3); % spike rate in spikes/s
spikeRate.rest = nSpikes.rest./(length(restEpochsIndex)*si_img*10^-3); 
spikeRate.mvm = nSpikes.mvm./(length(mvmEpochsIndex)*si_img*10^-3);

% movement rate score
mvmRateScore = (spikeRate.mvm-spikeRate.rest)./(spikeRate.mvm+spikeRate.rest);

% ratio n spikes in movement / n spikes in total
mvmTotRatio = nSpikes.mvm./spikeRate.tot; 
restTotRatio = nSpikes.rest./spikeRate.tot;

% Save results
fileName_rateMvmRest=strcat(fileName,'_rateMvmRest.mat');
save(fileName_rateMvmRest,'fileName', 'si_img','spikenums','nSpikes','spikeRate','nCells','-v7.3');
clear fileName_rateMvmRest;

%% Movement modulated cells - method1: z-score
% define movement start-on and movement start-off cells

exampleCell = 1; % set example cell to be plotted

baselineOffset = -2;

nBins = 50;

PSTHwindow = 3000; % period before and after movement onset to use (ms)

PSTHwindow = round(PSTHwindow/si_img); %converted to data points

mvmStartOnCells = zeros(size(spikenums.tot.orig,1),1);

mvmStartZscores = zeros(size(spikenums.tot.orig,1),nBins-1); %number of columns should be equal to nBins input argument of responseToEvent function called below

% Check z-score of example cell
figure;
[~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.tot.orig(exampleCell,:),mvmOnsetIndex,PSTHwindow,nBins,si_img,'post','plotting',baselineOffset); %#ok<ASGLU>
subplot(2,1,1)
title('Original data')

clear excitResponse inhibitResponse spikeZscore;


for cellCounter = 1:size(spikenums.tot.orig,1)
    [~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.tot.orig(cellCounter,:),mvmOnsetIndex,PSTHwindow,nBins,si_img,'post','noPlotting',baselineOffset);
    mvmStartOnCells(cellCounter)=excitResponse;
    mvmStartOffCells(cellCounter)=inhibResponse;
    mvmStartZscores(cellCounter,:) = spikeZscore;
end

mvmStartOnCells=find(mvmStartOnCells==true);
mvmStartOffCells=find(mvmStartOffCells==true);

clear excitResponse inhibitResponse spikeZscore;

    % reshuffling control
[spikenums.tot.shuffled] = circshiftRaster (spikenums.tot.orig,1000); %surrogate raster

% Check reshuffled z-score of example cell
figure;
[~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.tot.shuffled(exampleCell,:),mvmOnsetIndex,PSTHwindow,nBins,si_img,'post','plotting',baselineOffset);
subplot(2,1,2)
title('Surrogate data')

        % if the z-score of a cell is found significant in reshuffled spike count, the
        % cell is removed from the responsive group
for cellCounter = 1:size(spikenums.tot.orig,1)
    [~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.tot.shuffled(cellCounter,:),mvmOnsetIndex,PSTHwindow,nBins,si_img,'post','noPlotting',baselineOffset);
    if excitResponse == true && ismember(cellCounter,mvmStartOnCells)==1
        mvmStartOnCells(mvmStartOnCells==cellCounter)=[];
    end
    
    if inhibResponse == true && ismember(cellCounter,mvmStartOffCells)==1
        mvmStartOffCells(mvmStartOffCells==cellCounter)=[];
    end
end

clear excitResponse inhibResponse spikeZscore;

% define movement stop-on and movement stop-off cells
mvmStopOnCells = zeros(size(spikenums.tot.orig,1),1);

mvmStopZscores = zeros(size(spikenums.tot.orig,1),30); %number of columns should be equal to nBins input argument of responseToEvent function called below

for cellCounter = 1:size(spikenums.tot.orig,1)
    [~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.tot.orig(cellCounter,:),mvmOffsetIndex,PSTHwindow,nBins,si_img,'post','noPlotting');
    mvmStopOnCells(cellCounter)=excitResponse;
    mvmStopOffCells(cellCounter)=inhibResponse;
    mvmStopZscores(cellCounter,:) = spikeZscore;
end

mvmStopOnCells=find(mvmStopOnCells==true);
mvmStopOffCells=find(mvmStopOffCells==true);
mvmStopOffCells=mvmStopOffCells';

clear excitResponse inhibitResponse spikeZscore;

    % reshuffling control
%[spikenums.tot.shuffled] = circshiftRaster (spikenums.tot.orig,1000); %surrogate raster

        % if the z-score of a cell is found significant in reshuffled spike count, the
        % cell is removed from the responsive group
for cellCounter = 1:size(spikenums.tot.orig,1) 
    [~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.tot.shuffled(cellCounter,:),mvmOffsetIndex,PSTHwindow,nBins,si_img,'post','noPlotting');
    if excitResponse == true && ismember(cellCounter,mvmStopOnCells)==1
        mvmStopOnCells(mvmStopOnCells==cellCounter)=[];
    end
    
    if inhibResponse == true && ismember(cellCounter,mvmStopOffCells)==1
        mvmStopOffCells(mvmStopOffCells==cellCounter)=[];
    end
end

% group all responsive cells in struct array
%respCells.mvmModCellsZscore = unique([mvmStartOnCells; mvmStartOffCells; mvmStopOnCells; mvmStopOffCells']);
respCells.mvmStartOnCells = mvmStartOnCells;
respCells.mvmStartOffCells = mvmStartOffCells;
respCells.mvmStopOnCells = mvmStopOnCells;
respCells.mvmStopOffCells = mvmStopOffCells;

clear excitResponse inhibResponse spikeZscore mvmStartOnCells mvmStartOffCells mvmStopOnCells mvmStopOffCells;

% Save variables
if exist('fileName') == 0
    fileName = input('Type filename:','s');
end

fileName_mvmZscore=strcat(fileName,'_mvmZscore.mat');
save(fileName_mvmZscore,'fileName', 'spikenums','nSpikes','spikeRate','mvmStartZscores','mvmStopZscores','respCells','-v7.3');
clear fileName_mvmZscore;

%% Movement-modulated cells - method2: spike ratio and reshuffling
% spike rates in surrogate spike raster
if isfield(spikenums.tot,'shuffled')==0
    [spikenums.tot.shuffled] = circshiftRaster (spikenums.tot.orig,1000); %surrogate raster
end

nSpikes.tot_shuff = sum(spikenums.tot.shuffled,2); % n of spikes
nSpikes.rest_shuff = sum(spikenums.tot.shuffled(:,restEpochsIndex),2);
nSpikes.mvm_shuff = sum(spikenums.tot.shuffled(:,mvmEpochsIndex),2);

spikeRate.tot_shuff = nSpikes.tot_shuff./(size(spikenums.tot.shuffled,2)*si_img*10^-3); % spike rate in spikes/s
spikeRate.rest_shuff = nSpikes.rest_shuff./(length(restEpochsIndex)*si_img*10^-3); 
spikeRate.mvm_shuff = nSpikes.mvm_shuff./(length(mvmEpochsIndex)*si_img*10^-3);

% ratio n spikes in movement / n spikes in total (reshuffled spike raster)
mvmTotRatio_shuff = nSpikes.mvm_shuff./spikeRate.tot_shuff; 
restTotRatio_shuff = nSpikes.rest_shuff./spikeRate.tot_shuff;


% find thresholds to define movement and rest activated cells
mvmHighThresh = prctile(mvmTotRatio_shuff,99);
%restHighThresh = prctile(restTotRatio_shuff,95);

mvmLowThresh = prctile(mvmTotRatio_shuff,1);
%restLowThresh = prctile(restTotRatio_shuff,5);

[mvmOnCells] = find(mvmTotRatio>mvmHighThresh);
%[restOnCells] = find(restTotRatio>restHighThresh)

[mvmOffCells] = find(mvmTotRatio<mvmLowThresh);
%[restOffCells] = find(restTotRatio<restLowThresh)

% group all responsive cells in struct array
respCells.mvmOnCellsGrouped = unique([mvmOnCells; respCells.mvmStartOnCells]);
respCells.mvmModCellsRatio = unique([mvmOnCells; mvmOffCells]);
respCells.mvmModCellsGrouped = unique([respCells.mvmModCellsZscore; respCells.mvmModCellsRatio]);
respCells.mvmOnCells = mvmOnCells;
respCells.mvmOffCells = mvmOffCells;

clear mvmHighThresh mvmLowThresh mvmOnCells mvmOffCells

% Save variables
fileName_mvmRestMod=strcat(fileName,'_mvmRestMod.mat');
save(fileName_mvmRestMod,'fileName', 'spikenums','nSpikes','spikeRate','respCells','-v7.3');
clear fileName_mvmRestMod;

%% Definition of SCEs

SCEbinning = 10; % amount of data points to merge (default:4 --> 200ms) 
si_SCE = si_img * SCEbinning;

% Whole recording
clear rasterOptions; % options for plotting the raster plot
rasterOptions.cellColour = nCells.class2; % cells to highlight in binned spike raster
rasterOptions.epochs.epochs1 = 'SCE'; % highlight SCEs in binned spike raster
rasterOptions.epochs.epochs2 = mvmEpochsIndex./2; % highlight movement epochs in binned spike raster
[SCE.tot.nCellsSCE,SCE.tot.SCEindex,SCE.tot.SCE_count,nSpikes.SCE.tot,SCE.tot.SCEthresh,spikenums.tot.binned] = defSCEs(spikenums.tot.orig,si_SCE,1000,SCEbinning,rasterOptions); % define and plot SCEs
title('Spike raster, SCEs (green) and movement (grey)');


% Rest only
clear rasterOptions; % options for plotting the raster plot
rasterOptions.cellColour = nCells.class2; % cells to highlight in binned spike raster
rasterOptions.epochs.epochs1 = 'SCE'; % highlight SCEs in binned spike raster
[SCE.rest.nCellsSCE,SCE.rest.SCEindex,SCE.rest.SCE_count,nSpikes.SCE.rest,SCE.rest.SCEthresh,spikenums.rest.binned] = defSCEs(spikenums.tot.orig(:,restEpochsIndex),si_SCE,1000,SCEbinning,rasterOptions); % define and plot SCEs
title('Spike raster and SCEs (green) during rest');

% Cell participation in SCEs (rest only)
SCE.rest.cellsInSCE= zeros(size(spikenums.rest.binned,1), size(SCE.rest.SCEindex,2));
for cellCounter=1:size(spikenums.rest.binned,1)
    for SCEcounter=1:size(SCE.rest.SCEindex,2)
        if spikenums.rest.binned(cellCounter, SCE.rest.SCEindex (SCEcounter))==1
           SCE.rest.cellsInSCE(cellCounter,SCEcounter)=1;
        end
    end
end
SCE.rest.nSCEbyCell = sum(SCE.rest.cellsInSCE,2); % number of SCEs to which each cell participates
SCE.rest.fractSCEbyCell = SCE.rest.nSCEbyCell / length(SCE.rest.SCEindex); % proportion of SCEs to which each cell participates
figure;histogram(SCE.rest.fractSCEbyCell,20)

spikeRate.SCE.rest = nSpikes.SCE.rest ./(length(SCE.rest.SCEindex)*si_img*10^-3);

%% Define SCE-modulated cells (z-score)
% Only rest SCEs used

exampleCell = 1; % set example cell to be plotted

% define SCE-on/SCE-off cells
SCEonCells = zeros(size(spikenums.rest.binned,1),1);
SCEoffCells = zeros(size(spikenums.rest.binned,1),1);

SCEzScores = zeros(size(spikenums.rest.binned,1),21); %number of columns should be equal to nBins input argument of responseToEvent function called below

% Check z-score of example cell
figure;
[~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.rest.binned(exampleCell,:),SCE.rest.SCEindex,10,21,si_SCE,'during','plotting');
subplot(2,1,1)
title('Original data')


clear spikeZscore excitResponse inhibitResponse;

% Calculate z-scores and responses in original data
for cellCounter = nCells.tot
    [~,~,spikeZscore,signPosBins,~,excitResponse,inhibResponse] = responseToEvent (spikenums.rest.binned(cellCounter,:),SCE.rest.SCEindex,10,21,si_SCE,'during','noPlotting');
    SCEonCells(cellCounter)=excitResponse;
    SCEoffCells(cellCounter)=inhibResponse;
    SCEzScores(cellCounter,1:length(spikeZscore)) = spikeZscore;
end

SCEonCells=find(SCEonCells==true);
SCEoffCells=find(SCEoffCells==true);

clear spikeZscore excitResponse inhibitResponse;

%surrogate binned raster
[spikenums.rest.binned_shuffled] = circshiftRaster (spikenums.rest.binned,1000); %surrogate raster

% Check reshuffled z-score of example cell
figure;
[~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.rest.binned_shuffled(exampleCell,:),SCE.rest.SCEindex,10,21,si_SCE,'during','plotting');
subplot(2,1,2)
title('Surrogate data')


clear spikeZscore excitResponse inhibitResponse;

% Calculate z-scores and responses in surrogate data
for cellCounter = nCells.tot % if the z-score of a cell is found significant in reshuffled spike count, the
                                       % cell is removed from the responsive group
    [~,~,~,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.rest.binned_shuffled(cellCounter,:),SCE.rest.SCEindex,10,21,si_SCE,'during','noPlotting');
    if excitResponse == true && ismember(cellCounter,SCEonCells)==1
        SCEonCells(SCEonCells==cellCounter)=[];
    end
    
    if inhibResponse == true && ismember(cellCounter,SCEoffCells)==1
        SCEoffCells(SCEoffCells==cellCounter)=[];
    end
end

% group responsive cells in struct array
respCells.mvmPlusSCEonCells = intersect(respCells.mvmOnCellsGrouped,SCEonCells); %cells activated by both movement and SCEs
respCells.SCEonCells = SCEonCells;
respCells.SCEoffCells = SCEoffCells;

clear excitResponse inhibitResponse SCEcounter SCEonCells SCEoffCells;

% Save variables
fileName_SCE=strcat(fileName,'_SCE.mat');
save(fileName_SCE,'fileName', 'SCE','SCEbinning','respCells','spikenums','-v7.3');
clear fileName_SCE;

%% Definition of RUN-SCEs

SCEbinning = 1; % amount of data points to merge (default:4 --> 200ms) 
si_SCE = si_img * SCEbinning;

% Locomotion only
clear rasterOptions; % options for plotting the raster plot
rasterOptions.cellColour = [1]; % cells to highlight in binned spike raster
rasterOptions.epochs.epochs1 = 'SCE'; % highlight SCEs in binned spike raster
[SCE.mvm.nCellsSCE,SCE.mvm.SCEindex,SCE.mvm.SCE_count,nSpikes.SCE.mvm,SCE.mvm.SCEthresh,spikenums.mvm.binned] = defSCEs(spikenums.tot.orig(:,mvmEpochsIndex),si_SCE,1000,SCEbinning,rasterOptions); % define and plot SCEs
title('Spike raster and SCEs (green) during rest');

% Cell participation in SCEs (rest only)
SCE.mvm.cellsInSCE= zeros(size(spikenums.mvm.binned,1), size(SCE.mvm.SCEindex,2));
for cellCounter=1:size(spikenums.mvm.binned,1)
    for SCEcounter=1:size(SCE.mvm.SCEindex,2)
        if spikenums.mvm.binned(cellCounter, SCE.mvm.SCEindex (SCEcounter))==1
           SCE.mvm.cellsInSCE(cellCounter,SCEcounter)=1;
        end
    end
end
SCE.mvm.nSCEbyCell = sum(SCE.mvm.cellsInSCE,2); % number of SCEs to which each cell participates
SCE.mvm.fractSCEbyCell = SCE.mvm.nSCEbyCell / length(SCE.mvm.SCEindex); % proportion of SCEs to which each cell participates
figure;histogram(SCE.mvm.fractSCEbyCell,20)
ylabel('Number of cells')
xlabel('SCE participation rate / cell')

spikeRate.SCE.mvm = nSpikes.SCE.mvm ./(length(SCE.mvm.SCEindex)*si_img*10^-3);

%% Define RUN-SCE-modulated cells (z-score)
% Only rest SCEs used

exampleCell = 1; % set example cell to be plotted

% define SCE-on/SCE-off cells
runSCEonCells = zeros(size(spikenums.mvm.binned,1),1);
runSCEoffCells = zeros(size(spikenums.mvm.binned,1),1);

runSCEzScores = zeros(size(spikenums.mvm.binned,1),21); %number of columns should be equal to nBins input argument of responseToEvent function called below

% Check z-score of example cell
figure;
[~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.mvm.binned(exampleCell,:),SCE.mvm.SCEindex,10,21,si_SCE,'during','plotting');
subplot(2,1,1)
title('Original data')


clear spikeZscore excitResponse inhibitResponse;

% Calculate z-scores and responses in original data
for cellCounter = nCells.tot
    [~,~,spikeZscore,signPosBins,~,excitResponse,inhibResponse] = responseToEvent (spikenums.mvm.binned(cellCounter,:),SCE.mvm.SCEindex,10,21,si_SCE,'during','noPlotting');
    runSCEonCells(cellCounter)=excitResponse;
    runSCEoffCells(cellCounter)=inhibResponse;
    runSCEzScores(cellCounter,1:length(spikeZscore)) = spikeZscore;
end

runSCEonCells=find(runSCEonCells==true);
runSCEoffCells=find(runSCEoffCells==true);

clear spikeZscore excitResponse inhibitResponse;

%surrogate binned raster
[spikenums.mvm.binned_shuffled] = circshiftRaster (spikenums.mvm.binned,1000); %surrogate raster

% Check reshuffled z-score of example cell
figure;
[~,~,spikeZscore,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.mvm.binned_shuffled(exampleCell,:),SCE.mvm.SCEindex,10,21,si_SCE,'during','plotting');
subplot(2,1,2)
title('Surrogate data')


clear spikeZscore excitResponse inhibitResponse;

% Calculate z-scores and responses in surrogate data
for cellCounter = nCells.tot; % if the z-score of a cell is found significant in reshuffled spike count, the
                                       % cell is removed from the responsive group
    [~,~,~,~,~,excitResponse,inhibResponse] = responseToEvent (spikenums.mvm.binned_shuffled(cellCounter,:),SCE.mvm.SCEindex,10,21,si_SCE,'during','noPlotting');
    if excitResponse == true && ismember(cellCounter,runSCEonCells)==1
        runSCEonCells(runSCEonCells==cellCounter)=[];
    end
    
    if inhibResponse == true && ismember(cellCounter,runSCEoffCells)==1
        runSCEoffCells(runSCEoffCells==cellCounter)=[];
    end
end

% group responsive cells in struct array
respCells.mvmPlusSCEonCells = intersect(respCells.mvmOnCellsGrouped,SCEonCells); %cells activated by both movement and SCEs
respCells.runSCEonCells = runSCEonCells;
respCells.runSCEoffCells = runSCEoffCells;

clear excitResponse inhibitResponse SCEcounter SCEonCells SCEoffCells;

% Save variables
fileName_SCE=strcat(fileName,'_SCE.mat');
save(fileName_SCE,'fileName', 'SCE','SCEbinning','respCells','spikenums','-v7.3');
clear fileName_SCE;



