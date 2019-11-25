%--------------------------------------------------------------------------
% Spike extraction from 2-photon large-scale calcium imaging data
%%--------------------------------------------------------------------------

% >>> STRUCTURE >>>
%   -  rigid and non-rigid motion correction routine (NoRMCorre-master, https://github.com/flatironinstitute/NoRMCorre)
%   -  CaImAn pipeline:
%       1. constrained non-negative matrix factorization (segmentation and
%           trace extraction)
%       2. Markov chain Montecarlo for model-based spike inference
%   
% >>> REQUIRED THIRD-PARTY CODES >>>
%   -   NoRMCorre (https://github.com/flatironinstitute/NoRMCorre
%   -   CaImAn (https://github.com/flatironinstitute/CaImAn-MATLAB)
%   -   plot_staggered (Marco Bilucaglia)
%   -   sepblockfun (Matt Jacobson)
%
%--------------------------------------------------------------------------
% Yannick Bollman, Robin Dard, Marco Bocchio
%--------------------------------------------------------------------------
% Revisions
% - 2/8/19 Gaussian filter and detrending before trace selection
% 
%% LOAD DATA

%Load calcium movie
clear
gcp;

Y = concMovies;
Y = double(Y);      % convert to double precision 

% Plot mean of original movie
figure(1)
subplot(2,2,1)
imagesc(mean(Y,3))
title('Original movie')

%Load red channel image (if any)
sameNpixels = false; %same or different number of pixels between red and GCaMP channels?

isThereRedChannel = input('Is there a red channel? Y/N [Y]:','s');
if isempty(isThereRedChannel)
    isThereRedChannel = 'Y';
end
if isThereRedChannel=='Y';
    disp('select red channel image')
    [FileName2,PathName2,FilterIndex2]=uigetfile('*.tif');
    addpath(PathName2);
    redChannel = read_file(FileName2); % read the file
    redChSlice = input('Select optical slice to select from red channel: ','s');
    redChannel = imadjust(mat2gray(double(redChannel(:,:,str2num(redChSlice)))));     %adjust contrast and set double precision
    mapRedChannel = [0.3,0,0;0.4,0,0;0.5,0,0;0.6,0,0;0.7,0,0;0.8,0,0;0.9,0,0;1,0,0]; %colormap to visualise red channel;
    figure(2)
    imshow(redChannel,'Colormap',mapRedChannel);
end
    
    
    


%% PARAMETERS

refineCells = true; % Refine components (add/reject cells based on average matrix)
categCells =  1; % Accept/reject/categorize cells based on calcium trace
frameBlockAvg = 5; %number of frames to be averaged (see line 124)
nTrimmedPixels = 5; %number of pixels trimmed from the borders to crop movie after motion correction
fr = 9.62; %frame rate (in Hz)
si_img = 1/fr*10^3; %sampling interval (in ms)
[d1,d2,T] = size(Y);     


%% RIGID MOTION CORRECTION

%rigid motion correction
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',10,'max_shift',30,'us_fac',50,'init_batch',200,'iter',2);
tic; [M1,shifts1,template1] = normcorre(Y,options_rigid); toc

% Plot result
figure(1)
subplot(2,2,2)
imagesc(mean(M1,3))
title('After rigid motion correction')

%% NON-RIGID MOTION CORRECTION

%motion correction non-rigid
options_nonrigid = NoRMCorreSetParms('plot_flag',true,'d1',size(M1,1),'d2',size(M1,2),'grid_size',[20,20,1],'overlap_pre',[5,5,15],'mot_uf',[4,4,1],'bin_width',50,'max_shift',[5,5,5],'max_dev',[3,3,1],'us_fac',20,'init_batch',100,'boundary','copy','iter',1,'use_parallel',true,'buffer_width',150,'correct_bidir',false,'shifts_method','cubic');
tic; [M2,shifts2,template2] = normcorre(M1,options_nonrigid); toc
%clearvars M3

% Plot result
figure(1)
subplot(2,2,3)
imagesc(mean(M2,3))
title('After non-rigid motion correction')

%% ANISOTROPY CORRECTION (Arnaud Malvache)

Win = 1:500; % Temporal window used as reference
M3 = CorAnisotropy(M2, Win);
%figure(2)
imagesc(mean(M3,3))

% Plot mean of anisotropic filtered movie
figure(1)
subplot(2,2,4)
imagesc(mean(M3,3))
title('After anisotropic filtering')
  

%% RESIZING AND CROPPING


M3_cropped = cropMovie(M3,nTrimmedPixels); %crop movie to get rid of black borders caused by motion corrections (specify amount of pixels to be trimmed from each side)

if isThereRedChannel=='Y';
    redGreenPixelRatio = size(redChannel,1)/size(M3,1);    % Determine difference between number of pixels of calcium movie and of red channel
    if sameNpixels == false %resize red channel if nPixels different from GCaMP channel
    redChannel = imresize (redChannel, (1/redGreenPixelRatio)); %resize red channel image to match size of calcium image
    redChannelCropped = redChannel(1+nTrimmedPixels:end-nTrimmedPixels,1+nTrimmedPixels:end-nTrimmedPixels); %cropping of red channel image (same extent as calcium movie) 
    else
        redChannelCropped = redChannel;
    end
    imshow(redChannelCropped,'Colormap',mapRedChannel)
end

implay(mat2gray(M3_cropped));

[d1,d2,T] = size(M3_cropped);                                % dimensions of dataset
d = d1*d2;% total number of pixels



%% SAVE CORRECTED MOVIE

FileName_exist = exist ('fileName');
if FileName_exist == 0;
    prompt = 'Type filename: ';
    fileName = input(prompt);
end
    

fileName_CorrectedMovie=strcat(fileName,'_Corrected_Movie.mat');
save(fileName_CorrectedMovie,'fileName', 'M3_cropped','shifts2', 'redChannelCropped','fr','si_img','categCells','nTrimmedPixels','refineCells','-v7.3')

clearvars -except d d1 d2 fileName categCells fr fs_abf M3_cropped isThereRedChannel redChannel redChannelCropped refineCells si_abf si_img;

%% EXTRACTION PARAMETERS

K = 350;                                         % number of components to be found
tau = 4;                                        % std of gaussian kernel (half size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.75;                                  % merging threshold

options = CNMFSetParms(...   
    'fr',fr,...                                 %frame rate in Hz
    'decay_time',1,...                        %decay time in seconds (fs*decay time should be physical decay of transient)
    'ssub',1,...
    'tsub',2,...
    'd1',d1,'d2',d2,...                         % dimensionality of the FOV
    'p',p,...                                   % order of AR dynamics    
    'gSig',tau,...                              % half size of neuron
    'merge_thr',merge_thr,...                   % merging threshold  
    'deconv_method','constrained_foopsi',...
    'nb',1,...                                  % number of background components    
    'min_SNR',3,...                             % minimum SNR threshold
    'space_thresh',0.2,...                      % space correlation threshold
    'cnn_thr',0.2...                            % threshold for CNN classifier (not used)    
    );


if exist('fileName') == 0;
    fileName = input('Type filename:','s');
end

fileName_CaImAn_opts=strcat(fileName,'_CaImAn_opts.mat');
save(fileName_CaImAn_opts,'options', 'K', 'si_img');
clear fileName_CaImAn_opts


%% Data pre-processing

[P,Y] = preprocess_data(M3_cropped,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

% display centers of found components
Cn =  correlation_image(Y); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;


%% manually refine components (optional)
refineCells = true;  % flag for manual refinement
if refineCells
    [Ain,Cin,center]=manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end

%% update spatial components
Yr = reshape(Y,d,T);
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% classify components

%rval_space = classify_comp_corr(Y,A,C,b,f,options);
%ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test


%% event exceptionality

%fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
%ind_exc = (fitness < options.min_fitness);

%keep = ind_corr & ind_exc;

%A_keep = A(:,keep);
%C_keep = C(keep,:);

[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

Pm.p = 2;    % restore AR value
[A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);



%% EXTRACT TRACES
[A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F_new(A_or,C_or,b2,f2,P_or,options); % extract DF/F values (optional)
%add option
C_df=full(C_df);
[Coor,json_file] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints

%% CHECK TRACES
tracesToPlot = [1:size(C_df,1)];

plot_staggered(C_df(tracesToPlot,:));

clear tracesToPlot;

%% SMOOTH AND DETREND TRACES (optional)
tracesToFilter = [1:size(C_df,1)];

C_df_unfilt = C_df;
C_df(tracesToFilter) = detrend(C_df(tracesToFilter)); % linear detrend of traces
C_df(tracesToFilter) = smoothdata(C_df(tracesToFilter),'gaussian'); % smooth traces


%% MANUAL SELECTION ON TRACES 
%make trace selection
if refineCells == true
i=1;
%close all;

ContoursAll={};
%ContoursSoma={};
%ContoursPeri={};
ContoursRej={};
ContoursClass1={};
ContoursClass2={};

Tracecells=[];
Tracecellsdf=[];
%TraceNeuro=[];
%TraceNeurodf=[];
TraceNoise=[];
TraceNoisedf=[];
TraceClass1=[];
TraceClass1df=[];
TraceClass2=[];
TraceClass2df=[];

MemberShip=[];
BW4=zeros(d1,d2);
BW5=zeros(d1,d2);


selFigure=figure('units','normalized','outerposition',[0.25 0 0.5 1]);
hold on;

categ=0;
while i<=length(C_or(:,1))
[d1,d2] = size(Cn);
cont = medfilt1(Coor{i}')';
Cn2=Cn;
BW1=zeros(d1,d2);
BW2=zeros(d1,d2);
for j=1:size(cont,2)
BW1(cont(2,j),cont(1,j))=1;
 se1 = strel('line',2,90);
    se2 = strel('line',2,0);
    BWtemp=imdilate(BW1,[se1,se2]);
        BWtemp= imfill(BWtemp,'holes');
 BW2=bwperim(BWtemp);
end
subplot(6,2,1)
subplot(6,2,[1 2 3 4 5 6 7 8 9 10])
calciumImageSel = adapthisteq(Cn2); % calcium image to accept/reject cells

if isThereRedChannel=='Y';
compositeImage = imfuse(calciumImageSel,redChannelCropped,'ColorChannels','green-magenta'); %composite of red channel and calcium images
imagesc(compositeImage)
else
   imagesc(calciumImageSel)
end
hold on;

B = bwboundaries(BW2); %draw boundary of selected cell
visboundaries(B,'LineWidth',1,'Color','w');

if categ==28 %if 'go back' option was chosen at the last iteration
    if categHist(end)==30 %if last selected cell was in  category1 (uparrow choice, usually pyramidal cells)
    BW4=BW4-BW1; %remove last boundary added to category1 (affects only plotting)
    end
    
    if categHist(end)==31 %if last selected cell was in  category2 (downarrow choice, usually interneurons)
    BW5=BW5-BW1; %remove last boundary added to category2 (affects only plotting)
    end
   
end

B2=bwboundaries(BW4); %draw boundaries of cells in category1 (pyramidal cells)
h=visboundaries(B2,'LineWidth',1,'Color','b');
B3=bwboundaries(BW5); %draw boundaries of cells in category2 (interneurons)
visboundaries(B3,'LineWidth',1,'Color','r');


subplot(6,2,[11 12])
plot(C_df(i,:));
CellNum=num2str(i);
Ntot=num2str(size(C_df,1));
title(strcat('Cell #  ',CellNum, '/ ', Ntot))

k = waitforbuttonpress;
% 28 leftarrow
% 29 rightarrow
% 30 uparrow
% 31 downarrow
categ = double(get(gcf,'CurrentCharacter'));

if i==1
    categHist = categ;
end

if categ==30 %up arrow: set as pyramid
    Tracecells(end+1,:)=C_or(i,:);
    Tracecellsdf(end+1,:)=C_df(i,:);
    TraceClass1(end+1,:)=C_or(i,:);
    TraceClass1df(end+1,:)=C_df(i,:);
    ContoursAll{end+1}=cont;
    ContoursClass1{end+1}=cont;
    BW4=BW1+BW4;
end

if categ==31 %down arrow: set as interneuron
    Tracecells(end+1,:)=C_or(i,:);
    Tracecellsdf(end+1,:)=C_df(i,:);
    TraceClass2(end+1,:)=C_or(i,:);
    TraceClass2df(end+1,:)=C_df(i,:);
    ContoursAll{end+1}=cont;
    ContoursClass2{end+1}=cont;
    BW5=BW1+BW5;
end

if categ==29 %right arrow: reject cell
    TraceNoise(end+1,:)=C_or(i,:);
    TraceNoisedf(end+1,:)=C_df(i,:);
    ContoursRej{end+1}=cont;
end

if i>1
categHist = [categHist, categ]; %list of previous categorizations for detected cells
end

if categ==28 %if 'go back' option was chosen at the last iteration
    categHist=categHist(1:end-1);
    if categHist(end)==30 %if previous cell was in category1, remove all stored traces and contours from that cell
        Tracecells=Tracecells(1:end-1,:);
        Tracecellsdf=Tracecellsdf(1:end-1,:);
        TraceClass1=TraceClass1(1:end-1,:);
        TraceClass1df=TraceClass1df(1:end-1,:);
        ContoursAll=ContoursAll(:,1:end-1);
        ContoursClass1=ContoursClass1(:,1:end-1);         
    end
    
    if categHist(end)==31; %if previous cell was in category2, remove all stored traces and contours from that cell
        Tracecells=Tracecells(1:end-1,:);
        Tracecellsdf=Tracecellsdf(1:end-1,:);
        TraceClass2=TraceClass2(1:end-1,:);
        TraceClass2df=TraceClass2df(1:end-1,:);
        ContoursAll=ContoursAll(:,1:end-1);
        ContoursClass2=ContoursClass2(:,1:end-1);
    end
    
    if categHist(end)==29; %if previous cell was rejected, remove all stored traces and contours from that cell
        TraceNoise=TraceNoise(1:end-1,:);
        TraceNoisedf=TraceNoisedf(1:end-1,:);
        ContoursRej=ContoursRej(:,1:end-1);
    end             
        i=i-2;
    end
    
MemberShip(i)=categ;
i=i+1;
end
end 

MemberShip_filt=MemberShip;
MemberShip_filt(MemberShip_filt==2)=[];

ind_class1=transpose(find(MemberShip_filt==1));
ind_class2=transpose(find(MemberShip_filt==3));

C_df=Tracecellsdf;
%C_or=TraceNeuro;
%Coor=ContoursSoma;
%ContClass1=ContoursClass1;
%ContClass2=ContoursClass2;

Tracecells = [TraceClass2; TraceClass1]; %merge traces of the two populations (sparser cells first)
Tracecellsdf = [TraceClass2df; TraceClass1df]; %merge traces of the two populations (sparser cells first)


fileName_Traces=strcat(fileName,'_Traces.mat');
save(fileName_Traces,'Tracecells', 'Tracecellsdf', 'TraceClass1','TraceClass1df','TraceClass2','TraceClass2df','MemberShip','MemberShip_filt','ind_class1','ind_class2','-v7.3');

fileName_CellDetect=strcat(fileName,'_CellDetect.mat');
save(fileName_CellDetect,'d1','d2','ContoursAll','ContoursClass1','ContoursClass2','fr','si_img','-mat','-v7.3');

%% CONTOUR MAP FIGURES 
% Check if a second class of cells exists
oneCellClass= isempty(ContoursClass2);

%Build image (averaging full movie)
%Im=(imadjust(mat2gray(mean(Y,3))));
Im=(mat2gray(mean(M3_cropped,3)));

% Find position
pos.class1=zeros(length(ContoursClass1),2);
for i=1:length(ContoursClass1)
    PLS1=cell2mat(ContoursClass1(i));
    pos.class1(i,1)=median(PLS1(1,:));
    pos.class1(i,2)=median(PLS1(2,:));
end
pos.class1=round(pos.class1);

if oneCellClass==0 %if there are two cell classes
pos.class2=zeros(length(ContoursClass2),2);
for i=1:length(ContoursClass2)
    PLS2=cell2mat(ContoursClass2(i));
    pos.class2(i,1)=median(PLS2(1,:));
    pos.class2(i,2)=median(PLS2(2,:));
end
pos.class2=round(pos.class2);

pos.all = [pos.class2; pos.class1];

end

% ----------------------------------------------------
% Contours only on mean image
figure;
imagesc(Im);
%imagesc(mean(Y,3));
colormap(gray);
title('Contour map on mean image')
hold on;

for i=1:length(ContoursClass1)
cont = medfilt1(ContoursClass1{i}')';
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
B = bwboundaries(BW2,'noholes');
visboundaries(B,'Color','g','LineWidth', 0.2)
end

BW1=zeros(d1,d2);
BW2=zeros(d1,d2);

for i=1:length(ContoursClass2)
cont_class2 = medfilt1(ContoursClass2{i}')';
for j=1:size(cont_class2,2)
 BW1(cont_class2(2,j),cont_class2(1,j))=1;
 BW2=bwperim(BW1);
end
B = bwboundaries(BW2, 'noholes');
visboundaries(B,'Color','r','LineWidth', 0.2)
end
hold off;
savefig('Contour Map mean')



% ----------------------------------------------------
% Contours and cell numbers on mean image
figure;
imagesc(Im);
%imagesc(mean(Y,3));
colormap(gray);
title('Contour map on mean image')
hold on;

for i=1:length(ContoursClass1)
cont = medfilt1(ContoursClass1{i}')';
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
B = bwboundaries(BW2,'noholes');
visboundaries(B,'Color','g','LineWidth', 0.2)
end

BW1=zeros(d1,d2);
BW2=zeros(d1,d2);

for i=1:length(ContoursClass2)
cont_class2 = medfilt1(ContoursClass2{i}')';
for j=1:size(cont_class2,2)
 BW1(cont_class2(2,j),cont_class2(1,j))=1;
 BW2=bwperim(BW1);
end
B = bwboundaries(BW2, 'noholes');
visboundaries(B,'Color','r','LineWidth', 0.2)
end

% Add numbers
if oneCellClass==0 %if there are two cell classes
   
    for cellCounter=1:size(pos.class2,1)
        txt = num2str(cellCounter);
        text(pos.class2(cellCounter,1), pos.class2(cellCounter,2),txt,'Color','red','FontSize',12);
    end

    for cellCounter=1:size(pos.class1,1)
        txt = num2str(cellCounter+length(pos.class2));
        text(pos.class1(cellCounter,1), pos.class1(cellCounter,2),txt,'Color','green','FontSize',12);
    end
else
    for cellCounter=1:size(pos,1)
        txt = num2str(cellCounter);
        text(pos(cellCounter,1), pos(cellCounter,2),txt,'Color','green','FontSize',12);
    end

end
    
hold off;
savefig('Contour Map with numbers mean')
    


% ----------------------------------------------------
% Contours and numbers on correlation image
figure;
imagesc(calciumImageSel);
colormap(gray);
title('Contour map on correlation image')
hold on;

% Plot contours of class1
for i=1:length(ContoursClass1)
cont = medfilt1(ContoursClass1{i}')';
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
B = bwboundaries(BW2,'noholes');
visboundaries(B,'Color','g','LineWidth', 0.2)
end

BW1=zeros(d1,d2);
BW2=zeros(d1,d2);

if oneCellClass==0 %if there are two cell classes
    
    % Plot contours of class2
    for i=1:length(ContoursClass2)
    cont_class2 = medfilt1(ContoursClass2{i}')';
    for j=1:size(cont_class2,2)
     BW1(cont_class2(2,j),cont_class2(1,j))=1;
     BW2=bwperim(BW1);
    end
    B = bwboundaries(BW2, 'noholes');
    visboundaries(B,'Color','r','LineWidth', 0.2)
    end

% Add numbers      
    for cellCounter=1:size(pos.class2,1)
        txt = num2str(cellCounter);
        text(pos.class2(cellCounter,1), pos.class2(cellCounter,2),txt,'Color','red','FontSize',12);
    end

    for cellCounter=1:size(pos.class1,1)
        txt = num2str(cellCounter+length(pos.class2));
        text(pos.class1(cellCounter,1), pos.class1(cellCounter,2),txt,'Color','green','FontSize',12);
    end
else
    for cellCounter=1:size(pos,1)
        txt = num2str(cellCounter);
        text(pos(cellCounter,1), pos(cellCounter,2),txt,'Color','green','FontSize',12);
    end

end
    
hold off;

savefig('Contour Map with numbers correlation')

% ----------------------------------------------------
% Black & white contour map with filled contours
figure;

title('Contour map')

for i=1:length(ContoursAll)
cont = medfilt1(ContoursAll{i}')';
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
B = imfill(BW4,'holes');
imagesc(B)
map = [0 0 0;1 1 1];   % black & white color
colormap(map);
end

savefig('Contour Map filled BW')

%Colour-coded contour map
BW1=zeros(d1,d2);
BW2=zeros(d1,d2);

for i=1:length(ContoursClass1)
cont = medfilt1(ContoursClass1{i}')';
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
B=imfill(BW2,'holes');
map = [0 0 0;0 1 0];   % black and green color
colormap(map);
end

BW1=zeros(d1,d2);
BW2=zeros(d1,d2);

if oneCellClass==0 %if there are two cell classes
    for i=1:length(ContoursClass2)
    cont_class2 = medfilt1(ContoursClass2{i}')';
    for j=1:size(cont_class2,2)
     BW1(cont_class2(2,j),cont_class2(1,j))=1;
     BW2=bwperim(BW1);
    end
    Bi=imfill(BW2,'holes');
    map = [0 0 0;1 0 0];   % black and red color
    colormap(map);
    end

end

figure
title('Contour map filled')
ContMap=imfuse(B,Bi,'ColorChannels',[2 1 0]);
imagesc(ContMap)
savefig('Contour Map filled colours')

fileName_pos=strcat(fileName,'_pos.mat');
save(fileName_pos, 'pos')
clear fileName_pos;

   
    %% EXTRACT SPIKES (Class1)
       
    tic;
    disp(['MCMC on test cell running...']);
    
    %MCMC Class1
    [SAMPLES, defparams.class1] = cont_ca_sampler(TraceClass1df(20,:));    %% MCMC example 
    plot_continuous_samples(SAMPLES,TraceClass1df(20,:));

    disp(['MCMC on all cells running...']);

    
    att={};
    parfor i = 1:length(TraceClass1df(:,1))
    try
    SAMPLES = cont_ca_sampler(TraceClass1df(i,:));    %% MCMC
    att{i}=SAMPLES.params.spiketimes_;
    catch
        att{i}=0;
    end
    end

    %make Class1 spikes
    d1=length(TraceClass1df(:,1));
    d2=2*length(TraceClass1df(1,:));
    spikenumsClass1=zeros(d1,d2);
    for i = 1:length(TraceClass1df(:,1))
        ac=att{i};
        if ac>0
        acc = 0.5;
        result = round(ac/acc)*acc;
        result=result*2;
        for j=1:length(result)
            spikenumsClass1(i,result(j))=1;
        end
        end
    end
    
    attClass1 = att;
    
    
    disp(['MCMC done']);
toc;
    
%% EXTRACT SPIKES (Class2)
       
    tic;
    disp(['MCMC on test cell running...']);
    
    %MCMC Class2
    [SAMPLES, defparams.class2] = cont_ca_sampler(TraceClass2df(1,:));    %% MCMC example 
    plot_continuous_samples(SAMPLES,TraceClass2df(1,:));

    disp(['MCMC on all cells running...']);

    
    att={};
    parfor i = 1:length(TraceClass2df(:,1))
    try
    SAMPLES = cont_ca_sampler(TraceClass2df(i,:));    %% MCMC
    att{i}=SAMPLES.params.spiketimes_;
    catch
        att{i}=0;
    end
    end

    %make Class2 spikes
    d1=length(TraceClass2df(:,1));
    d2=2*length(TraceClass2df(1,:));
    spikenumsClass2=zeros(d1,d2);
    for i = 1:length(TraceClass2df(:,1))
        ac=att{i};
        if ac>0
        acc = 0.5;
        result = round(ac/acc)*acc;
        result=result*2;
        for j=1:length(result)
            spikenumsClass2(i,result(j))=1;
        end
        end
    end
    
   attClass2 = att;
   
    disp(['MCMC done']);
toc;

    spikenums = [spikenumsClass2; spikenumsClass1]; %merge spike raster of the two populations (sparser cells first)

   
    fileName_MCMC=strcat(fileName,'_MCMC');
    save(fileName_MCMC, 'spikenums', 'attClass1', 'spikenumsClass1', 'attClass2', 'spikenumsClass2', 'pos')
    
      
    
%% EXTRACT SPIKES (ALL CELLS)

tic;
    disp(['MCMC on test cell running...']);

%MCMC
[SAMPLES, defparams.all] = cont_ca_sampler(Tracecells(1,:));    %% MCMC example 
plot_continuous_samples(SAMPLES,Tracecells(1,:));

    disp(['MCMC on all cells running...']);


att={};
parfor i = 1:length(Tracecells(:,1))
try
SAMPLES = cont_ca_sampler(Tracecells(i,:));    %% MCMC
att{i}=SAMPLES.params.spiketimes_;
catch
    att{i}=0;
end
end

%make spikes
d1=length(Tracecells(:,1));
d2=2*length(Tracecells(1,:));
spikenums=zeros(d1,d2);
for i = 1:length(Tracecells(:,1))
    ac=att{i};
    if ac>0
    acc = 0.5;
    result = round(ac/acc)*acc;
    result=result*2;
    for j=1:length(result)
        spikenums(i,result(j))=1;
    end
    end
end

disp(['MCMC done']);
toc;

%% VERIFY DETECTION
plotRasterTraces (spikenums,Tracecellsdf,'bin');

%% SAVE RESULTS
nCells.tot = 1:size(spikenums,1);
nCells.class1 = size(spikenumsClass2,1)+1:size(spikenumsClass1,1);
nCells.class2 = 1:size(spikenumsClass2,1);
clear grp
grp = [zeros(1,nCells.class2),ones(1,nCells.class1)];

fileName_MCMC=strcat(fileName,'_MCMC.mat');
save(fileName_MCMC, 'att', 'spikenums','pos','defparams')
clear fileName_MCMC;

fileName_preprocSummary=strcat(fileName, '_preprocSummary.mat');
save(fileName_preprocSummary, 'Tracecells', 'Tracecellsdf', 'att', 'fileName', 'pos', 'spikenums', 'nCells', 'grp')
clear fileName_preprocSummary;




