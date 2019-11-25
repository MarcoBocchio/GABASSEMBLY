function XCells=NoCamp_Overlap(Cells,Nx,Ny)

XCells=Cells;
while true
    AllCells=sum(XCells,3);
    % Find first overlapped area
    [a,b]=ind2sub([Nx,Ny],find(AllCells>1,1));
    if isempty(a)
        break
    end
    % Find overlapped cells
    NumXCell=find(XCells(a,b,:)==1);
    
    tmp=sum(XCells(:,:,NumXCell),3);
    if sum(tmp(:)>1)>sum(tmp(:)==1)/3
        XCells(:,:,NumXCell(1))=tmp>0;
        XCells(:,:,NumXCell(2:end))=[];
    else
        for j=1:length(NumXCell)
            XCells(:,:,NumXCell(j))=XCells(:,:,NumXCell(j))-double(tmp>1);
        end
    end
end