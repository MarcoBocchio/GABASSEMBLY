% Continuous movement correction based on cross-correlation barycenter
% in a selected window
function ImgCor=CorMovWin2(Img,Win)

%Display mean image
h=figure;
imagesc(mean(Img(:,:,Win),3)) %Coming soon: automatic box
Box=ginput(2);
close(h)

%Define data in box
WinX=floor(min(Box(:,2))):floor(max(Box(:,2)));
WinY=floor(min(Box(:,1))):floor(max(Box(:,1)));
[Nx,Ny,Nz]=size(Img);
ImgBox=Img(WinX,WinY,:);

%Define reference and its FT
Ref=mean(ImgBox(:,:,Win),3);
TFRef=conj(fft2(Ref));

%Cross-correlation and shift detection
Shift=zeros(Nz,2);
for i=1:Nz
    Cor=(fftshift(ifft2(fft2(ImgBox(:,:,i)).*TFRef)));
    Shift(i,:)=bary(Norm01(Cor).^10);
end

%Shift image
dy=Shift(:,1)-length(WinY)/2;
dx=Shift(:,2)-length(WinX)/2;
parfor i=1:Nz
    ImgCor(:,:,i)=interp2((1:Ny)',1:Nx,Img(:,:,i),((1:Ny)+dy(i))',(1:Nx)+dx(i));
end