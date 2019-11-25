function [Fringe,Step]=NoCamp_ReadPos2(data)

%Use a threshold from data distribution
data0=GaussBlur1d(Norm01(data),length(data)/100,1);
th=graythresh(data0);
dataTh=data0>th;
tmp=diff(dataTh);
Fringe=find(abs(tmp)==1);
if Fringe(1)~=1
    Fringe=[1;Fringe];
end
Step=find(diff(Fringe)>1e4);