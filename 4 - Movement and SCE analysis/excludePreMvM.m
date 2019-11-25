for mvmIndexCounter = 1:length(mvmOnsetIndex)
preMvmIndexTemp = mvmOnsetIndex(mvmIndexCounter)-20:mvmOnsetIndex(mvmIndexCounter)-1;
if mvmIndexCounter == 1
preMvmIndex = preMvmIndexTemp;
else
preMvmIndex = [preMvmIndex preMvmIndexTemp];
end
end

restEpochsTrimmedLogic = restEpochsLogic;
restEpochsTrimmedLogic(preMvmIndex)=[];
restEpochsTrimmedIndex = find (restEpochsTrimmedLogic == 1);