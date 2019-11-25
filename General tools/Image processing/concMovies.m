function Z = loadMovie;
%%
% >>> OPERATION >>>
% Concatenation of short TIFF calcium movies into a long one

%%
% >>> USE >>>
%Select first TIFF file in the folder. The remaining files will be
%automatically found and concatenated

%% 
% Marco Bocchio, updated 27/6/19

clear

[FileName,PathName,FilterIndex]=uigetfile('*.tif');
addpath(PathName);

X = read_file(FileName);

numFramesSingle = size(X,3);

files = dir(PathName);
fileNames = {files.name};

fileSizes = {files.bytes};
fileSizes = (cell2mat(fileSizes))./1000000; %size of files in MB

tifFiles = (fileNames(fileSizes>20));

Y=zeros(size(X,1),size(X,2),numFramesSingle*length(tifFiles)+numFramesSingle-1);

tic;

for i=1:length(tifFiles)
   tifFile=cell2mat(tifFiles(:,i));
      Y = read_file(tifFile);
  
   if i==1
       Z = Y;
   else   
   Z = cat(3,Z,Y);
   end
   
end
  toc;
  
[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels
   


end

