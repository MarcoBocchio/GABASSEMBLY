% load movie information
[File0, Path0] = uigetfile('*.tif');
info=imfinfo([Path0,File0]);
Nz=length(info);
Ny=info(1).Width;
Nx=info(1).Height;
%A=zeros(Nx,Ny,Nz*NTif);

%% Load movie and first average

movie = zeros (Nx, Ny, Nz);
for i = 1: Nz
    movie (:, : , i) = imread (File0, 'TIFF', 'INDEX', i);
end

Ref= mean (movie,3);
figure
imagesc(Ref)
colormap (gray)

savefig ('Ref1.fig')

%% Correct Anisotropy

movie2=CorAnisotropy(movie, 1:Nz);
Ref2=mean(movie2,3);
figure
imagesc(Ref2)
colormap(gray)

savefig ('Ref2.fig')

%% Correct Movements: Arnaud's method

movie3=CorMovWin2(movie2, 1:Nz);
%movie3=movie2; %if already corrected for movements
Ref3= mean(movie3,3);
figure
imagesc(Ref3)
colormap (gray)

savefig ('Ref3.fig')


%% write MotCorreFinal tif movie

for i=1:Nz
    imwrite(uint16(movie3(:,:,i)),strcat(File0,'_Corre.tif'),'TIFF','writemode', 'append');
end

%% PCA-ICA

movie3=double(movie3);
movie3(isnan(movie3))=0;
movie3=double(movie3);
IC=NoCamp_PCA_ICA(movie3);

fil1D=IC.ic_XFilter';
TrIC=IC.ic_TimeCourse';
NCamp=size(fil1D,1);

%% COURSE 3

% remove the one that are not skewed enough, removing things (may be removed to analyze pups)
SkSp = skewness(fil1D'); %Spaial skewness
SkT = skewness(TrIC') ; % temporal skewness

%% Select GCamp traces
disp('Select traces')
SkSpMin = 1; % 1 for RACE or GCAMP6; 5 for place cells; the smaller the value the more we keep
CompOKtmp = find(abs(SkSp)>SkSpMin & SkT>1);
button = zeros(1,length(CompOKtmp));
figure
%Nx=200;
%Ny=200;
for i = 1:length(CompOKtmp)
    subplot(1,5,1:4)
    plot(TrIC(CompOKtmp(i),:))
    subplot(1,5,5)
    imagesc(reshape(fil1D(CompOKtmp(i),:),Nx,Ny))
    title([num2str(i),'/',num2str(length(CompOKtmp))])
    [~,~,button(i)] = ginput(1);
end

CompOK = CompOKtmp(button==1);
CompPerisom = CompOKtmp(button==32);

%% Find cells from covariant pixels
disp('Find cells and traces')
Cells = [];
NeuronSize=10;
PxSize =2;
for i = CompOK
    if SkSp(i)<-1 && SkT(i)<-1
        fil1D(i,:)=-fil1D(i,:);
    end
    CellsNew = NoCamp_FindCells(fil1D(i,:),Nx,Ny,NeuronSize,PxSize);
    Cells = cat(3,Cells,CellsNew);
end
clear CellsNew
NpxCells=squeeze(sum(sum(Cells)));
% Remove small and big areas
NeurArea=pi*(NeuronSize/PxSize/2)^2;
CellOK=NpxCells>NeurArea/2 & NpxCells<NeurArea*2;

%% Remove or merge overlapping areas
XCells=NoCamp_Overlap(Cells(:,:,CellOK),Nx,Ny);
NCell=size(XCells,3);

%% Close cells
se = strel('disk',5);
XCellsCl=zeros(Nx,Ny,NCell);
for i=1:NCell
    XCellsCl(:,:,i)=imclose(XCells(:,:,i),se);
end
%Remove empty cells
EmptyCell = [];
for i = 1:NCell
    if sum(sum(XCellsCl(:,:,i))) == 0
        EmptyCell = [EmptyCell i];
    end
end
XCellsCl(:,:,EmptyCell) = [];
NCell=size(XCellsCl,3);
figure
imagesc(sum(XCellsCl,3))
colormap gray

savefig ('Cell_Map.fig')

%% Calculate cells' GCAMP signal
Nz=14000;
Img1D=reshape(movie3,[Nx*Ny Nz]);
Tr0=zeros(NCell,Nz);
for i=1:NCell
    tmp=XCellsCl(:,:,i)>0;
    Tr0(i,:)=mean(Img1D(tmp,:));
end
clear tmp
Tr1=zeros(NCell,Nz);
for i=1:NCell
    Tr1(i,:)=Tr0(i,:)/median(Tr0(i,:));
end

%PLOT TRACES
figure
for i = 1 : NCell
    plot(Tr1 (i,:) + i-1)
    hold on
end

savefig ('Traces.fig')

save('Arnaud_part1','-v7.3')

