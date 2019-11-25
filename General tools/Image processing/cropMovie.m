function [croppedMovie] = cropMovie(inputMovie,croppedPixels)

nPixels = size(inputMovie,1);
nFrames = size(inputMovie,3);
nCroppedMoviePixels = nPixels - (croppedPixels*2);
croppedMovie = zeros(nCroppedMoviePixels,nCroppedMoviePixels,nFrames);

for i = 1:size (inputMovie,3);
    croppedMovie(:,:,i) = inputMovie(1+croppedPixels:end-croppedPixels,1+croppedPixels:end-croppedPixels,i);
end

end

