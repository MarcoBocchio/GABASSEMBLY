SkSp = skewness(Fil1d');
SkT = skewness(TrIC');

%% Select GCamp traces
disp('Select traces')
SkSpMin = 2; % 1 for RACE or GCAMP6; 5 for place cells
CompOKtmp = find(abs(SkSp)>SkSpMin & SkT>1);
button = zeros(1,length(CompOKtmp));
figure
Nx=200;
Ny=200;
for i = 1:length(CompOKtmp)
    subplot(1,5,1:4)
    plot(TrIC(CompOKtmp(i),:))
    subplot(1,5,5)
    imagesc(reshape(Fil1d(CompOKtmp(i),:),Nx,Ny))
    title([num2str(i),'/',num2str(length(CompOKtmp))])
    [~,~,button(i)] = ginput(1);
end

CompOK = CompOKtmp(button==1);
CompPerisom = CompOKtmp(button==32);

NeuronSize=10;
PxSize=2;

%% Find cells from covariant pixels
disp('Find cells and traces')
Cells = [];
for i = CompOK
    if SkSp(i)<-1 && SkT(i)<-1
        Fil1d(i,:)=-Fil1d(i,:);
    end
    CellsNew = NoCamp_FindCells(Fil1d(i,:),Nx,Ny,NeuronSize,PxSize);
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

Nz = 4000;

%% Calculate cells' GCAMP signal
Img1D=reshape(Movie3,[Nx*Ny Nz]);
Tr0=zeros(NCell,Nz);
for i=1:NCell
    tmp=XCellsCl(:,:,i)>0;
    Tr0(i,:)=mean(Img1D(tmp,:));
end
%if Smooth == 1
 %   Tr0 = Tr0 + circshift(Tr0,1,2); %Smooth movement
%end
clear tmp
Tr1=zeros(NCell,Nz);
for i=1:NCell
    Tr1(i,:)=Tr0(i,:)/median(Tr0(i,:));
end

%% Remove identical cells
[Tr1,XCellsCl] = NoCamp_MergeSameCellsFunc(Tr1,XCellsCl);
NCell = size(Tr1,1);

%% PhotoBleaching
Tr1b=zeros(NCell,Nz);
ws = warning('off','all');
for i=1:NCell
    p0=polyfit(1:Nz,Tr1(i,:),3);
    Tr1b(i,:)=Tr1(i,:)./polyval(p0,1:Nz);
end
warning(ws)

%% Display traces
figure
for i=1:NCell
    plot(MovT,Tr1b(i,:)+i-1)
    hold on
end

%% Save temporal data
cd(PathSave)
save([PathSave,'\MovT.mat'],'MovT')
save([PathSave,'\WinRest.mat'],'WinRest')
save([PathSave,'\Cells.mat'],'XCellsCl')
save([PathSave,'\Tr1b.mat'],'Tr1b')
save([PathSave,'\Ephy.mat'],'Ephy')
save([PathSave,'\FOV.mat'],'FOV')

%% Split running and rest
Sr = 2e4;
Nz = 1000;
NFilms = 2;

%% Load abf file
fnmabf = fullfile(pathname2,filenameAbf);
data = abfload(fnmabf);
clear fnmabf filenameAbf pathname2

tE = (1:length(data))/Sr; %Measurement time
Ephy = downsample(data(:,ChEphy),10);

%% Divide movies

MovT1 = NoCamp_MovieTime_Gaps(data1(:,ChTrig),Sr,Nz/NFilms);

%% Calculate LapXPT
LapXPT = interp1(tE,data(:,ChLapXP),MovT,'pchip');

%% Find run part and speed
%!!! New data
[Fringe,Step] = NoCamp_ReadPos2(data(:,ChSpeed));
NRun = length(Step)-1; RunLim = zeros(NRun,2); RestLim = zeros(NRun,2);
for i = 1:NRun
    RunLim(i,:) = [Fringe(Step(i)+1) Fringe(Step(i+1))];
    RestLim(i,:) = [Fringe(Step(i)) Fringe(Step(i)+1)];
end

%% Concatenate running sessions (wider window : -50/+50 frames)
Win = [];
for i = 1:length(RunLim)
    a = find(MovT>RunLim(i,1)/Sr,1);
    b = find(MovT<RunLim(i,2)/Sr,1,'last');
    Win = [Win max(a-50,1):min(b+50,Nz)];
end
Win = unique(Win);

%% Concatenate resting sessions (narrower window : +50/-50 frames)
WinRest = [];
for i = 1:length(RestLim)
    a = find(MovT>RestLim(i,1)/Sr,1);
    b = find(MovT<RestLim(i,2)/Sr,1,'last');
    WinRest = [WinRest a+50:b-50];
end
WinRest = unique(WinRest);


%% Distance and speed
disp('Spatial Domain')
NFringe=100;

while true
    % In Electrophy time
    DistE = interp1(Fringe/Sr,(1:length(Fringe))/NFringe/2,(1:max(Fringe))/Sr,'pchip'); %In laps
    DistE(isnan(DistE)) = max(DistE);
    %add last points (extrapolation) with small increase (to avoid
    %interpolation bug)
    DistE = [DistE max(DistE)+(1:length(tE)-length(DistE))*1e-10];
    SpeedE = diff(DistE)*LTreadmill*Sr; %In cm/s
    % In Movie time
    D = interp1(tE,DistE,MovT,'pchip');
    Speed = interp1(tE(2:end),SpeedE,MovT);
    Dist = min(D):0.001:max(D);
    LapXP = interp1(D,LapXPT,Dist);
    LapXP = [LapXP zeros(1,1e3-rem(length(LapXP),1e3))];
    figure(10)
    imagesc(reshape(-LapXP,1e3,floor(length(LapXP)/1e3))')
    title(['NFringe', num2str(NFringe)]);
    choice = questdlg('Are you happy ?','LapXP');
    switch choice
        case 'Yes'
            break
        case 'No'
            NFringe=str2double(inputdlg({'Fringes in 1 Lap'},'Adjust LapXP',1,{int2str(NFringe)}));
        case 'Cancel'
            break
    end
end
saveas(gcf,[PathSave,'\lapxp.fig'],'fig');saveas(gcf,[PathSave,'\lapxp.png'],'png');

%% Conversion in the spatial domain
TrSp=zeros(NCell,length(Dist));
for i=1:NCell 
    if isnan(Tr1b(i,1))
        Tr1b(i,:)=0;
    end
    TrSp(i,:)=interp1(D,Tr1b(i,:),Dist,'pchip');
    TrSp(i,:)=Norm01(TrSp(i,:));
end
SpeedSp=interp1(D,Speed,Dist);

%% Save Spatial Data

save([PathSave,'\TrSp.mat'],'TrSp')
save([PathSave,'\Speed.mat'],'Speed')
save([PathSave,'\SpeedSp.mat'],'SpeedSp')
save([PathSave,'\Laps.mat'],'D')
save([PathSave,'\Distance.mat'],'Dist')