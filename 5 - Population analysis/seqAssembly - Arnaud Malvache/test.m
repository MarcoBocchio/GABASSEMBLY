for i=1:NCamp
subplot (1,5,1:4)
plot (TrIC(i,:))
subplot(1,5,5)
imagesc(reshape (fil1D(i,:),192,191))
ginput (1);
end