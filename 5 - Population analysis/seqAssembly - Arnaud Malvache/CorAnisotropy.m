function D = CorAnisotropy(Img0,Win)

[Nx,Ny,Nz] = size(Img0);
%Mean over steady frames
A = mean(Img0(:,:,Win),3);
%Padding to avoid wrapping
B = [flipud([fliplr(A) A fliplr(A)]);[fliplr(A) A fliplr(A)];flipud([fliplr(A) A fliplr(A)])];
%Low Frequency Pass
C = GaussBlur(B,Nx/20);
%Apply correction
D = single(zeros(Nx,Ny,Nz));
for i = 1:Nz
    D(:,:,i) = Img0(:,:,i)./sqrt(C(Nx+1:2*Nx,Ny+1:2*Ny));
end