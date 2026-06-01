function [anovaMat]=anovaOptoAndDelay(allDataTestsOnly,allDataTestsOnlyDelay,anovaMat,anovaMatDelay)
%%% d' stats
testFullTrial=NaN(18,2); 
testFullTrial(1:6,1:2)=allDataTestsOnly{2,39}';
testFullTrial(7:12,1:2)=allDataTestsOnly{2,40}';
testFullTrial(13:18,1:2)=allDataTestsOnly{2,41}';

conditions=5;
linIdx=1:length(firstIdxAnimal):(length(firstIdxAnimal)*conditions);
testFullTrialD=NaN(length(firstIdxAnimal)*conditions,2); 
valIdx=39;
for ii=1:length(linIdx)
    testFullTrialD(linIdx(ii):linIdx(ii)+length(firstIdxAnimal)-1,1)=allDataTestsOnlyDelay{2,valIdx}(1,:);
    testFullTrialD(linIdx(ii):linIdx(ii)+length(firstIdxAnimal)-1,2)=allDataTestsOnlyDelay{2,valIdx}(2,:);
    valIdx=valIdx+2;
end
subsetDelay=downsample(testFullTrialD,2);
testFullTrial(19:24,1:2)=subsetDelay(1:6,:);
% 6 animals per condition for the original data, so now downsample the
% delay data to 6 samples and add to the matrix

[~,~,statsMGBbyAnDprime]=anova2(testFullTrial,6);
c1An = multcompare(statsMGBbyAnDprime);
tbl1An = array2table(c1An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2An = multcompare(statsMGBbyAnDprime,"Estimate","row");
tbl2AnD = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

%%% hit stats
hitStats=NaN(18,2); 
hitStats(1:6,1)=allDataTestsOnly{6,39}(1,1:6)';
hitStats(1:6,2)=allDataTestsOnly{6,39}(3,1:6)';
hitStats(7:12,1)=allDataTestsOnly{6,40}(1,1:6)';
hitStats(7:12,2)=allDataTestsOnly{6,40}(3,1:6)';
hitStats(13:18,1)=allDataTestsOnly{6,41}(1,1:6)';
hitStats(13:18,2)=allDataTestsOnly{6,41}(3,1:6)';

hitStatsD=NaN(length(firstIdxAnimal)*conditions,2); 
valIdx=18;
for tt=1:conditions
    for yy=1:length(firstIdxAnimal)
        delayTemp=allDataTestsOnlyDelay{firstIdxAnimal(yy),valIdx};
        reinfTemp=allDataTestsOnlyDelay{firstIdxAnimal(yy),16};
        delayTemp=cell2mat(delayTemp);reinfTemp=cell2mat(reinfTemp);
        reinfTemp(isnan(delayTemp))=nan;
        delayVal(yy,:)=delayTemp;
        reinfVal(yy,:)=reinfTemp;
    end
    valIdx=valIdx+2;
    if tt==1
        delayValAnimal=delayVal;
        reinfValAnimal=reinfVal;
    else
        delayValAnimal=vertcat(delayValAnimal,delayVal); % each row is an animal, every 3 rows is a condition
        reinfValAnimal=vertcat(reinfValAnimal,reinfVal);
    end
end
delayValAvgAnimal=nanmean(delayValAnimal,2);
reinfValAvgAnimal=nanmean(reinfValAnimal,2);
hitStatsD=horzcat(reinfValAvgAnimal,delayValAvgAnimal);
subsetDelay=downsample(hitStatsD,2);
hitStats(19:24,1:2)=subsetDelay(1:6,:);

[~,~,statsMGBbyAnHit]=anova2(hitStats,6);
c1An = multcompare(statsMGBbyAnHit);
tbl1An = array2table(c1An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2An = multcompare(statsMGBbyAnHit,"Estimate","row");
tbl2AnHit = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);


%%% fa stats
FAStats=NaN(18,2); 
FAStats(1:6,1)=allDataTestsOnly{6,39}(2,1:6)';
FAStats(1:6,2)=allDataTestsOnly{6,39}(4,1:6)';
FAStats(7:12,1)=allDataTestsOnly{6,40}(2,1:6)';
FAStats(7:12,2)=allDataTestsOnly{6,40}(4,1:6)';
FAStats(13:18,1)=allDataTestsOnly{6,41}(2,1:6)';
FAStats(13:18,2)=allDataTestsOnly{6,41}(4,1:6)';

FAStatsD=NaN(length(firstIdxAnimal)*conditions,2); 
valIdx=19;
for tt=1:conditions
    for yy=1:length(firstIdxAnimal)
        delayTemp=allDataTestsOnlyDelay{firstIdxAnimal(yy),valIdx};
        reinfTemp=allDataTestsOnlyDelay{firstIdxAnimal(yy),17};
        delayTemp=cell2mat(delayTemp);reinfTemp=cell2mat(reinfTemp);
        reinfTemp(isnan(delayTemp))=nan;
        delayVal(yy,:)=delayTemp;
        reinfVal(yy,:)=reinfTemp;
    end
    valIdx=valIdx+2;
    if tt==1
        delayValAnimal=delayVal;
        reinfValAnimal=reinfVal;
    else
        delayValAnimal=vertcat(delayValAnimal,delayVal); % each row is an animal, every 3 rows is a condition
        reinfValAnimal=vertcat(reinfValAnimal,reinfVal);
    end
end
delayValAvgAnimal=nanmean(delayValAnimal,2);
reinfValAvgAnimal=nanmean(reinfValAnimal,2);
FAStatsD=horzcat(reinfValAvgAnimal,delayValAvgAnimal);
subsetDelay=downsample(FAStatsD,2);
FAStats(19:24,1:2)=subsetDelay(1:6,:);

[~,~,statsMGBbyAnFA]=anova2(FAStats,6);
c1An = multcompare(statsMGBbyAnFA);
tbl1An = array2table(c1An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2An = multcompare(statsMGBbyAnFA,"Estimate","row");
tbl2AnFA = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
