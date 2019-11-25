function [imagingTrimmed, mvmTrimmed, fs_abf, si_abf,timeStart,timeEnd] = loadMvm(spikenums,imagingCh, mvmCh, timeStart,timeEnd);

% >>> OPERATION >>>
% Load abf file with imaging-on, treadmill (and LFP if applicable) data.
% Option: choose specific interval of abf to isolate (e.g. to analyse only
% one movie of the series)

% >>> INPUT VARIABLES >>>
% NAME             TYPE, DEFAULT        DESCRIPTION
% spikenums        double               spike raster matrix
% imagingCh        scalar               imaging channel # on abf file
% mvmCh            scalar               locomotion channel # on abf file
% timeStart        scalar               start time of temporal window to isolate a specific part of the abf (optional)         
% timeEnd          scalar               start time of the temporal window (optional)

% >>> REQUIRED THIRD-PARTY CODES >>>
% - abfload (F. Collman)

% 
% Marco Bocchio, 5/7/19

%% load ABF
[fileName,pathName] = uigetfile('*.abf');
fullFileName = fullfile(pathName,fileName);
[abfTrace,h,si]=abfload (fullFileName);
si_abf = si.si/10^3; %sampling interval (in ms)
fs_abf = 1/si_abf; %sampling rate (in kHz)

%% isolate specific time window
if nargin > 3
    pntStart = round(timeStart*60/si_abf*10^3); %start point    WRONG: CORRECT
    pntEnd = round(timeEnd*60/si_abf*10^3); %end point          WRONG:CORRECT
    abfTrace = abfTrace(pntStart:pntEnd,:);
end

%% load channels
imaging=abfTrace(:,imagingCh);
mvm = abfTrace(:,mvmCh);

%% isolate 'imaging on' epochs
[imagingOffIndex] = consecAboveThresh(imaging,0.02,1000);
imagingTrimmed=imaging;
imagingTrimmed([imagingOffIndex])=[];
mvmTrimmed=mvm;
mvmTrimmed([imagingOffIndex])=[];

%% plot
close all;
%spike raster
subplot(3,1,1)
plotSpikeRaster(logical(spikenums),'PlotType','vertline');
xlim([0 length(spikenums)])

%movement
subplot(3,1,2)
plot(mvmTrimmed);
xlim([0 length(mvmTrimmed)])

%imaging on
subplot(3,1,3)
plot(imagingTrimmed);
xlim([0 length(imagingTrimmed)])

end






