function SpGrL=WaveletLine2D(A,f0,sigma)  

[nx,ny]=size(A);
if rem(nx,2)==0
    kx=-nx/2+0.5:nx/2-0.5;
else
    kx=-(nx-1)/2:(nx-1)/2;
end
if rem(ny,2)==0
    ky=-ny/2+0.5:ny/2-0.5;
else
    ky=-(ny-1)/2:(ny-1)/2;
end
[Ky,Kx]=meshgrid(ky,kx);
Krho=sqrt(Kx.^2+Ky.^2);
% Fourier Gate (ring filter)
G=exp(-(Krho-f0).^2/sigma.^2)/sqrt(sigma); %Morlet 2D

TFA=fftshift(fft2(A));
% Fourier Transform
SpGrL=ifft2(ifftshift(TFA.*G));