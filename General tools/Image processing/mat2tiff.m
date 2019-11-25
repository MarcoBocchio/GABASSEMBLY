function mat2tiff(fileName, Y, frameStart,frameEnd)

tic;
%write movie to tiff

if nargin == 2;
    frameStart = 1;
    frameEnd = size(Y,3);
end

SS = fileName;
prefix='MotC_';
for frameCounter=frameStart:frameEnd
    %disp(num2str(frameCounter));
    imwrite(uint16(Y(:,:,frameCounter)), strcat(prefix,SS,'.tif'),'TIFF','writemode', 'append');
end

toc;


end