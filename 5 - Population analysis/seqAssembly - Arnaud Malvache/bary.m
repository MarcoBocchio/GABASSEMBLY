function b=bary(A)

[Nx,Ny]=size(A);
x=1:Nx;y=1:Ny;
[X,Y]=meshgrid(y,x);
bx=sum(sum(A.*X))/sum(A(:));
by=sum(sum(A.*Y))/sum(A(:));
b=[bx,by];