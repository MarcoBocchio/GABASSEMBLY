function [index, logicArray] = consecAboveThresh (inputVector, amplThreshold, spanThreshold)

%% >>> OPERATION >>>
% Finds n consecutive values above a certain threshold

%% 
% Marco Bocchio, 4/7/19

if size(inputVector,2) == 1;
    inputVector = inputVector';
end
    

aboveThreshold = (inputVector > amplThreshold);  %where above threshold
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);  
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth > spanThreshold;   
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
index = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));
logicArray = zeros(length(inputVector),1);
logicArray(index)=1;

end
