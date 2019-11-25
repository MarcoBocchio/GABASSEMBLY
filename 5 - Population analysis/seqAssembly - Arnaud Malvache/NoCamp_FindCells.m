%% Detect Cells from PCA-ICA for a given component

function Cells=NoCamp_FindCells(A0,Nx,Ny,NeuronSize,PxSize)

% Get 2D filter
A=reshape(A0,Nx,Ny);

% 2D Wavelet
f0=Nx/(NeuronSize/PxSize);
B=real(WaveletLine2D(A,f0,f0/2));

% Find extrema
[a,b,c,d]=extrema2(B);
Cell=b(a>max(B(:))/5);%3 for GCAMP5 arbitrary threshold
[dx,dy]=ind2sub([Nx,Ny],Cell);
Border=dx<7 | dx>Nx-7 | dy<7 | dy>Ny-7;
dx(Border)=[];
dy(Border)=[];
% figure(1)
% hold off
% imagesc(A)
% hold on
% plot(dy,dx,'+k')

% Find cells from each extrema
Cells=zeros(Nx,Ny,length(dx));
se = strel('disk',1);
L=bwlabel(imerode(B>max(B(:))/10,se));
for i=1:length(dx)
    Cells(:,:,i)=imdilate(L==L(dx(i),dy(i)),se);
end
% figure(2)
% imagesc(sum(Cells,3))
% colormap gray