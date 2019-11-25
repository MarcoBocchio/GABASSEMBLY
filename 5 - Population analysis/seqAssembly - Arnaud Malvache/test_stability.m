for i = 1:NClOK;
    
    eventsassemblie =find(PCl (i,:)==1);
    nbrCells = RCl(i,eventsassemblie);
   
    meannbrcells(i)= mean(nbrCells);
    prop(i)= meannbrcells(i)/size (C0{i},2);
     
end 

mean(prop)