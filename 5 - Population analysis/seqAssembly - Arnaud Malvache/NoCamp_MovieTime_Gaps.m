function MovT = NoCamp_MovieTime_Gaps(data,Sr,Nz)

%% old data (summer 2014)
% UP=[1;find(data>-1)];
% Step=find(diff(UP)>1000);
% for i=1:NFilms
%     if i+1>length(Step)
%         MovLim(i,:)=[UP(Step(i)+1) UP(end)];
%     else
%         MovLim(i,:)=[UP(Step(i)+1) UP(Step(i+1))];
%     end
%     dttmp=diff(MovLim(i,:))/Sr/Nz; %Duration of a movie frame
%     MovT(i,:)=MovLim(i,1)/Sr+(dttmp/2:dttmp:dttmp*Nz-dttmp/2); %Time (s) for each movie
% end

%% new data
% Fr = find(diff(data)>0.01);
% TFr = median(diff(Fr))/Sr;
% MovT = Fr(1)/Sr + (TFr/2:TFr:TFr*Nz-TFr/2); %Time (s) for all movies

%% new data with gaps
% Find trigs
Down = find(diff(data)<-0.01);
UP = find(diff(data)>0.01);
Steps_d = diff(Down);
Down(1) = []; %remove one point to align with diff
Steps_u = diff(UP);
Step_size = median(Steps_d);

Gaps_start = Down(Steps_d>3*Step_size);
Gaps_start = [Down(1);Gaps_start];
Gaps_end = UP(Steps_u>3*Step_size);
Gaps_end = [Gaps_end;UP(end)];
display(['Number of detected gaps: ',int2str(length(Gaps_start))])

MovT = [];
for i = 1:length(Gaps_start)
    Dt = (Gaps_end(i)-Gaps_start(i));
    MovT = [MovT Gaps_start(i):Dt/(Nz-1):Gaps_end(i)];
end
MovT = MovT/Sr;