function [Y_adj] = imadjust_movie(Y)

Y_adj = zeros(size(Y,1),size(Y,2),size(Y,3));

for i=1:size(Y,3)
Y_adj(:,:,i)=imadjust(Y(:,:,i));
   
end

