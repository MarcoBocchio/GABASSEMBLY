%clear all
%close all
clear
Nz=size(raw_traces,2);
NCell=size(raw_traces,1);
Tr=zeros(NCell,Nz);
for i=1:NCell
    Tr(i,:)=raw_traces(i,:)/median(raw_traces(i,:));
end

%load('Speed');

%%
[NCell,Nt]=size(Tr);
dt=200;
[M4,PCout] = step1_PCA(GaussBlur1d(Tr,Nt/20,2),10);
P=PCout.pc_TimeCourse';
%%
SigRef=diff(P(1,:));
% Calculate Shifted trace
Cor1d=zeros(1,NCell);
Shift1d=zeros(1,NCell);
Trsh=zeros(NCell,Nt);
for i=1:NCell
        tmp=covnorm(Tr(i,:),SigRef,dt);
        Cor1d(i)=max(tmp);
        Shift1d(i)=find(tmp==max(tmp))-dt-1;
        Trsh(i,:)=circshift(Tr(i,:)',-Shift1d(i))';
end

DC=find(Cor1d>graythresh(Cor1d));
NDC=length(DC);
TrDC=Tr(DC,:);
for i=1:NDC
    TrDC(i,:)=Norm01(GaussBlur1d(TrDC(i,:),Nt/2,2));
    TrDC(i,:)=Norm01(TrDC(i,:));
    TrDC(i,:)=TrDC(i,:)-median(TrDC(i,:));
end
ShiftDC=Shift1d(DC);
[~,xDel]=sort(ShiftDC);


Assembly=sum(Trsh(DC,:));

%Sp=downsample(Speed,10);
Sp=downsample(Speed,1);
SpBlur=GaussBlur(Sp,Nt/10);
t=(1:Nt)/10;
% t=MovT;

figure
subplot(4,1,1)
imagesc(t(1:floor(Nt/4)),1:NDC,TrDC(xDel,1:floor(Nt/4)),[0 0.5])
hold on
plot(t(1:floor(Nt/4)),NDC-SpBlur(1:floor(Nt/4))/2,'g')
subplot(4,1,2)
imagesc(t(floor(Nt/4):floor(Nt/2)),1:NDC,TrDC(xDel,floor(Nt/4):floor(Nt/2)),[0 0.5])
hold on
plot(t(floor(Nt/4):floor(Nt/2)),NDC-SpBlur(floor(Nt/4):floor(Nt/2))/2,'g')
subplot(4,1,3)
imagesc(t(floor(Nt/2):floor(3*Nt/4)),1:NDC,TrDC(xDel,floor(Nt/2):floor(3*Nt/4)),[0 0.5])
hold on
plot(t(floor(Nt/2):floor(3*Nt/4)),NDC-SpBlur(floor(Nt/2):floor(3*Nt/4))/2,'g')
subplot(4,1,4)
imagesc(t(floor(3*Nt/4):Nt),1:NDC,TrDC(xDel,floor(3*Nt/4):Nt),[0 0.5])
colormap hot
hold on
plot(t(floor(3*Nt/4):Nt),NDC-SpBlur(floor(3*Nt/4):Nt)/2,'g')
xlabel('Time (in s)')

%break
%% Savings
save('DC.mat','DC')
save('DCShift.mat','ShiftDC')
save('TrDC.mat','TrDC')