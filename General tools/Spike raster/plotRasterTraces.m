function plotRasterTraces (spikeRaster,traces,bin);

%% Downsample spike raster to match traces

N_cells = size(traces,1);
time = 1:size(traces,2);

switch bin
    case 'bin'
       spikenums_mod = logical(sepblockfun(spikeRaster, [1,2], 'max'));
        
    case 'noBin'
        spikenums_mod = logical(spikeRaster);
end

%imagesc(spikenums_mod)

Norm = zeros(size(spikenums_mod));



%% Plot the graph

for cellCounter=1:(N_cells)
    Norm(cellCounter,:) = traces(cellCounter,:)/max(traces(cellCounter,:));
end


 % old method
         %for cellCounter = 1:N_cells
            %plot(Norm(cellCounter,:) + 1*(cellCounter-1),'k')
            %plot(find(spikenums_mod(cellCounter,:)), 1*(cellCounter-1), '.','MarkerSize',5,'MarkerEdgeColor','r')  
          %end
   
    
  figure
  hold on;
     
   for cellCounter = 1:N_cells
        spikePos = time(spikenums_mod(cellCounter, :));
        plot(Norm(cellCounter,:) + cellCounter,'k')
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [cellCounter-0.1 cellCounter+0.5], 'r'); %else in black
        end
   end 
   
   set(gca,'children',flipud(get(gca,'children')))
   ylim([0 N_cells+1])
   xlim([0 length(time)])
   
   
end