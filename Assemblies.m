%--------------------------------------------------------------------------
% Pipeline to analyse assembly activities in large groups of imaged cells
%--------------------------------------------------------------------------

% >>> STRUCTURE >>>
%   ASSEMBLY DETECTION
%   -   seqAssembly set of codes for k-means clustering of SCEs for assembly detection (Malvache et al., 2016)
%   -   SGCassembly: similarity graph clustering (Avitan et al., 2017, Moelner et al., 2018)
%   -   ICAssembly: PCA/ICA method (Lopes Dos Santos et al., 2013)
%   -
% >>> REQUIRED THIRD-PARTY CODES >>>
%   -   toolbox (Vitor Lopes-Dos Santos, vtlsantos@gmail.com)
%   -   seqAssembly (Arnaud Malvache)
%   -   SGCAssembly (similarity graph clustering method, Jan Moelter, part of the neural assembly detection
%       package: https://github.com/GoodhillLab/neural-assembly-detection)

%--------------------------------------------------------------------------
% Marco Bocchio, updated 6/8/2019
%--------------------------------------------------------------------------
% Updates and revisions
% 06.08.19  fixed bugs with definition of single and multiple assembly cells 


%% Exclude cells from assembly detection

cellsToExclude = [];

spikenums.rest.assemblies = spikenums.rest.binned;
spikenums.rest.assemblies(cellsToExclude,:) = [];
spikenums.rest.assemblies_shuffled=spikenums.rest.binned_shuffled;
spikenums.rest.assemblies_shuffled(cellsToExclude,:)=[];

nCells.assemblies = size(spikenums.rest.assemblies,1);


%% SCEkm_assembly

% run SCEkm_assembly code
SCEkm_assembly;

% find all cells forming assemblies
SCEkm_results.cellsInAssemblies = sum(CellCl,1);
SCEkm_results.cellsInAssemblies = find(SCEkm_results.cellsInAssemblies == 1);
SCEkm_results.propCellsInAssemblies = length(SCEkm_results.cellsInAssemblies)./nCells.tot;

% save remaining results in struct array
SCEkm_results.NCl = NCl;clear NCl;
SCEkm_results.NClOK = NClOK;clear NClOK;
SCEkm_results.cellsInSCE = cellsInSCE;clear cellsInSCE;
SCEkm_results.x1 = x1;clear x1;
SCEkm_results.RList = RList;clear RList;
SCEkm_results.MSort = MSort;clear MSort;
SCEkm_results.CellCl = CellCl;clear CellCl;
SCEkm_results.assemblyMembers = C0;clear C0;
SCEkm_results.nAssemblies = length(SCEkm_results.assemblyMembers);
SCEkm_results.propCellsInAssemblies = length(SCEkm_results.cellsInAssemblies)./nCells.tot;

% assembly overlap (cells shared in multiple assemblies)
if SCEkm_results.nAssemblies > 1
    for assemblyCounter = 1:SCEkm_results.nAssemblies
        if assemblyCounter == 1
        concCellAssemblies = SCEkm_results.assemblyMembers{1,1};
        end
        if assemblyCounter > 1
        concCellAssemblies = horzcat(concCellAssemblies,SCEkm_results.assemblyMembers{1,assemblyCounter});
        end
    end

        %pre-allocation
    SCEkm_results.singleAssemblyCells = zeros(length(SCEkm_results.cellsInAssemblies),1);
    SCEkm_results.multiAssemblyCells = zeros(length(SCEkm_results.cellsInAssemblies),1);

        % find single assembly and multiple assemblies cells
    for cellAssemblyCounter = 1:length(SCEkm_results.cellsInAssemblies);
        cellInAssembly_idx = find(concCellAssemblies == SCEkm_results.cellsInAssemblies(cellAssemblyCounter));
        if length(cellInAssembly_idx) == 1
            SCEkm_results.singleAssemblyCells(cellAssemblyCounter,1) = SCEkm_results.cellsInAssemblies(cellAssemblyCounter);
        else
            SCEkm_results.multiAssemblyCells(cellAssemblyCounter,1) = SCEkm_results.cellsInAssemblies(cellAssemblyCounter);
        end
    end
          % remove zero elements
    SCEkm_results.singleAssemblyCells(SCEkm_results.singleAssemblyCells==0) = [];
    SCEkm_results.multiAssemblyCells(SCEkm_results.multiAssemblyCells==0) = [];

        % assembly overlap score (n multi assembly cells / n tot cells forming
        % assemblies)
    SCEkm_results.assemblyOverlap = length(SCEkm_results.multiAssemblyCells)./length(SCEkm_results.cellsInAssemblies);
else
    SCEkm_results.singleAssemblyCells = SCEkm_results.nAssemblies;
    SCEkm_results.multiAssemblyCells = [];
    SCEkm_results.assemblyOverlap = [];

end

% save results
fileName_SCEkm_assembly=strcat(fileName,'_SCEkm_assembly.mat');
save(fileName_SCEkm_assembly,'SCEkm_results','-v7.3');
clear fileName_SCEkm_assembly;

%% SGC_assembly

% rearrange and save arrays to run SGC code
activity_raster = spikenums.rest.binned';
activity_raster_peaks = SCE.rest.SCEindex;
activity_raster_peak_threshold = SCE.rest.SCEthresh;
save('_ACTIVITY-RASTER.mat','activity_raster','activity_raster_peaks','activity_raster_peak_threshold');

% run SGC code
[ SGC_results ] = SGC_ASSEMBLY_DETECTION('_ACTIVITY-RASTER.mat')

% number of cells forming assemblies
SGC_results.nAssemblies = length(SGC_results.assemblies);
SGC_results.cellsInAssemblies = SGC_results.assemblies{1,1};
if SGC_results.cellsInAssemblies > 1
    for assemblyCounter = 2:SGC_results.nAssemblies
        SGC_results.cellsInAssemblies = horzcat(SGC_results.cellsInAssemblies, SGC_results.assemblies{assemblyCounter,1});
    end
end
SGC_results.cellsInAssemblies = unique(SGC_results.cellsInAssemblies');
SGC_results.propCellsInAssemblies = length(SGC_results.cellsInAssemblies)./nCells.assemblies;

% assembly overlap (cells shared in multiple assemblies)
    %concatenate cells forming assemblies in one vector
for assemblyCounter = 1:SGC_results.nAssemblies
    if assemblyCounter == 1
    concCellAssemblies = SGC_results.assemblies{1,1};
    end
    if assemblyCounter > 1
    concCellAssemblies = horzcat(concCellAssemblies,SGC_results.assemblies{assemblyCounter,1});
    end
end

    %pre-allocation
SGC_results.singleAssemblyCells = zeros(length(SGC_results.cellsInAssemblies),1);
SGC_results.multiAssemblyCells = zeros(length(SGC_results.cellsInAssemblies),1);

    % find single assembly and multiple assemblies cells
for cellAssemblyCounter = 1:length(SGC_results.cellsInAssemblies);
    cellInAssembly_idx = find(concCellAssemblies == SGC_results.cellsInAssemblies(cellAssemblyCounter));
    if length(cellInAssembly_idx) == 1
        SGC_results.singleAssemblyCells(cellAssemblyCounter,1) = SGC_results.cellsInAssemblies(cellAssemblyCounter);
    else
        SGC_results.multiAssemblyCells(cellAssemblyCounter,1) = SGC_results.cellsInAssemblies(cellAssemblyCounter);
    end
end
      % remove zero elements
    SGC_results.singleAssemblyCells(SGC_results.singleAssemblyCells==0) = [];
    SGC_results.multiAssemblyCells(SGC_results.multiAssemblyCells==0) = [];

    % assembly overlap score (n multi assembly cells / n tot cells forming
    % assemblies)
SGC_results.assemblyOverlap = length(SGC_results.multiAssemblyCells)./length(SGC_results.cellsInAssemblies);
    
clear uniqueCellAssemblies_idx concCellAssemblies

% find data points in which assemblies are active
for assemblyCounter = 1:SGC_results.nAssemblies
    for SCEcounter = 1:length(SGC_results.assembly_pattern_detection.assemblyIActivityPatterns{assemblyCounter,1})
        SGC_results.assemblyActivities{assemblyCounter,1}(1,SCEcounter) = SCE.rest.SCEindex(SGC_results.assembly_pattern_detection.assemblyIActivityPatterns{assemblyCounter,1}(SCEcounter));
    end
end

SGC_results.assemblyActivLogic = zeros (SGC_results.nAssemblies, length(spikenums.rest.binned));
for assemblyCounter = 1:SGC_results.nAssemblies
    for SCEcounter = 1:length(SGC_results.assembly_pattern_detection.assemblyIActivityPatterns{assemblyCounter,1})
        SCEindex = (SGC_results.assemblyActivities{assemblyCounter,1}(SCEcounter));
    SGC_results.assemblyActivLogic(assemblyCounter,SCEindex) = 1;
    end
end
clear SCEindex SCEcounter

% Cells/assemblies activity probabilities
for assemblyCounter = 1:SGC_results.nAssemblies
    for cellCounter=1:nCells.tot
        assembly_cell = ([logical(spikenums.rest.binned(cellCounter,:)); SGC_results.assemblyActivLogic(assemblyCounter,:)]);
        assembly_cell_coact = (sum(assembly_cell,1));
        SGC_results.P_assembly_cell(cellCounter,assemblyCounter) = length(find(sum(assembly_cell,1)==2))./length(SGC_results.assemblyActivities{assemblyCounter,1}); % probability that a given cell fires when a given assembly is active
        SGC_results.P_cell_Assembly(cellCounter,assemblyCounter) = length(find(sum(assembly_cell,1)==2))./length(att{1,cellCounter}); % probability that a given assembly is active when a given cell fires
    end
end

clear cellCounter assemblyCounter

% save results
fileName_SGC_assembly=strcat(fileName,'_SGC_assembly.mat');
save(fileName_SGC_assembly,'SGC_results','-v7.3');

clear fileName_SGC_input fileName_SGC_assembly;

%% ICA_assembly

%spikenums.rest.binned = spikenums.tot.binned(:,restEpochsIndex);
%[spikenums.rest.shuffled] = circshiftRaster (spikenums.rest.binned,1000,'noPlotting');

 
% Plot correlation matrix
ICA_results.correlationmat = corr(spikenums.rest.binned');
figure
imagesc(ICA_results.correlationmat)

    %Options
clear opts;
opts.threshold.permutations_percentile = 99;
opts.threshold.number_of_permutations = 1000;
opts.Patterns.number_of_iterations = 500;
opts.threshold.method = 'circularshift';
opts.Patterns.method = 'ICA';

    % Assembly patterns
ICA_results.Patterns = assembly_patterns(spikenums.rest.binned,opts);
ICA_results.nAssemblies = size(ICA_results.Patterns,2);

   
    %Assembly activities
ICA_results.Activities = assembly_activity(ICA_results.Patterns,spikenums.rest.binned);

if exist('assemblyMembers')==1
    clear assemblyMembers
end

ICA_results.cellsInAssemblies = [];

figure;
% Cell participation and plot of assembly patterns
for assemblyCounter = 1:ICA_results.nAssemblies
   % Find cells that are part of an assembly (>1.5 standard deviation of
        % weights)
    highestWeightSign = abs(max(ICA_results.Patterns(:,assemblyCounter)))-abs(min(ICA_results.Patterns(:,assemblyCounter)));
    if highestWeightSign > 0
        assemblyMembers_temp = find(zscore(ICA_results.Patterns(:,assemblyCounter))>2);
    else
        assemblyMembers_temp = find(zscore(ICA_results.Patterns(:,assemblyCounter))<-2);
    end
    ICA_results.assemblyMembers{assemblyCounter} = assemblyMembers_temp;
    ICA_results.cellsInAssemblies = vertcat(ICA_results.cellsInAssemblies, assemblyMembers_temp);
    
        % Plot assembly weights
    subplot(ICA_results.nAssemblies,1,assemblyCounter)
    stem(ICA_results.Patterns(:,assemblyCounter))
end

ICA_results.cellsInAssemblies = unique(ICA_results.cellsInAssemblies);
ICA_results.propCellsInAssemblies = length(ICA_results.cellsInAssemblies)./nCells.assemblies;

         
% assembly overlap (cells shared in multiple assemblies)
    %concatenate cells forming assemblies in one vector
for assemblyCounter = 1:ICA_results.nAssemblies
    if assemblyCounter == 1
    concCellAssemblies = ICA_results.assemblyMembers{1,1};
    end
    if assemblyCounter > 1
    concCellAssemblies = vertcat(concCellAssemblies,ICA_results.assemblyMembers{1,assemblyCounter});
    end
end

    %pre-allocation
ICA_results.singleAssemblyCells = zeros(length(ICA_results.cellsInAssemblies),1);
ICA_results.multiAssemblyCells = zeros(length(ICA_results.cellsInAssemblies),1);

    % find single assembly and multiple assemblies cells
for cellAssemblyCounter = 1:length(ICA_results.cellsInAssemblies);
    cellInAssembly_idx = find(concCellAssemblies == ICA_results.cellsInAssemblies(cellAssemblyCounter));
    if length(cellInAssembly_idx) == 1
        ICA_results.singleAssemblyCells(cellAssemblyCounter,1) = ICA_results.cellsInAssemblies(cellAssemblyCounter);
    else
        ICA_results.multiAssemblyCells(cellAssemblyCounter,1) = ICA_results.cellsInAssemblies(cellAssemblyCounter);
    end
end
      % remove zero elements
ICA_results.singleAssemblyCells(ICA_results.singleAssemblyCells==0) = [];
ICA_results.multiAssemblyCells(ICA_results.multiAssemblyCells==0) = [];

    % assembly overlap score (n multi assembly cells / n tot cells forming
    % assemblies)
ICA_results.assemblyOverlap = length(ICA_results.multiAssemblyCells)./length(ICA_results.cellsInAssemblies);
    
clear uniqueCellAssemblies_idx concCellAssemblies

% Plot graphs
figure
subplot(211)
imagesc(spikenums.rest.binned)
xlim([0 length(spikenums.rest.binned)])
subplot(212)
plot(ICA_results.Activities')
xlim([0 length(spikenums.rest.binned)])
legend

% Save results
fileName_ICA_assembly=strcat(fileName,'_ICA_results.mat');
save(fileName_ICA_assembly,'fileName', 'ICA_results','-v7.3');
clear fileName_ICA_assembly;

%% Cross-correlation of single cell activities to assembly activities (ICA)

cellToPlot = 1;
maxLag = 20;

[ICA_xcorr_results.posCorrCells,ICA_xcorr_results.negCorrCells,ICA_xcorr_results.meanCorrCells,ICA_xcorr_results.corr_cellAssembly,ICA_xcorr_results.meanCorr_cellAssembly] = corrCellToAssembly(spikenums.rest.binned_all,spikenums.rest.assemblies_shuffled,ICA_results.Activities,maxLag,si_SCE,cellToPlot);

ICA_xcorr_results.nCorrAssemblies = sum(ICA_xcorr_results.posCorrCells,1);

ICA_xcorr_results.highRespCells = intersect(respCells.mvmPlusSCEonCells,ICA_xcorr_results.meanCorrCells);

fileName_ICA_assembly_xcorr=strcat(fileName,'_ICA_assembly_xcorr.mat');
save(fileName_ICA_assembly_xcorr,'fileName', 'ICA_xcorr_results','-v7.3');
clear fileName_ICA_assembly_xcorr;

