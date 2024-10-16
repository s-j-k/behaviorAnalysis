function [allDataTestsOnly,mgbTempTestsOnly,tempTestsOnly]=diffPcOptoSess(allDataTestsOnly,mgbTempTestsOnly,tempTestsOnly,optocolor)

wwFig=figure(100);
clear jj
for jj=2:size(mgbTempTestsOnly,1)
    allDataTestsOnly{jj,39}=mgbTempTestsOnly{jj,28}-mgbTempTestsOnly{jj,27};
    allDataTestsOnly{jj,40}=mgbTempTestsOnly{jj,30}-mgbTempTestsOnly{jj,29};
    allDataTestsOnly{jj,41}=mgbTempTestsOnly{jj,32}-mgbTempTestsOnly{jj,31};
end
diffFull=NaN; diffTone=NaN; diffChoice=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    diffFull=cat(1,diffFull,allDataTestsOnly{jj,39});
    diffTone=cat(1,diffTone,allDataTestsOnly{jj,40});
    diffChoice=cat(1,diffChoice,allDataTestsOnly{jj,41});
end
diffChoiceTrunc=diffChoice(~isnan(diffChoice));
diffChoiceTrunc=diffChoiceTrunc(1:length(diffFull(~isnan(diffFull))));
catMGBDiff=[diffFull(~isnan(diffFull)) diffTone(~isnan(diffTone)) diffChoiceTrunc];
subplot(1,2,1); hold on;
qqq=bar([nanmean(diffFull) nanmean(diffTone) nanmean(diffChoice)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffFull,1),1), ...
    diffFull,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffTone,1),1), ...
    diffTone,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffChoice,1),1), ...
    diffChoice,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title(['MGB Inactivation']);
xticklabels({'Full Trial','Tone','Choice'});
xlabel('Condition');

clear jj
for jj=2:size(tempTestsOnly,1) %IC
    allDataTestsOnly{jj,42}=tempTestsOnly{jj,34}-tempTestsOnly{jj,33};
    allDataTestsOnly{jj,43}=tempTestsOnly{jj,36}-tempTestsOnly{jj,35};
    allDataTestsOnly{jj,44}=tempTestsOnly{jj,38}-tempTestsOnly{jj,37};
end
diffFull=NaN; diffTone=NaN; diffChoice=NaN;
for jj=2:size(tempTestsOnly,1)
    diffFull=cat(1,diffFull,allDataTestsOnly{jj,42});
    diffTone=cat(1,diffTone,allDataTestsOnly{jj,43});
    diffChoice=cat(1,diffChoice,allDataTestsOnly{jj,44});
end
diffChoiceTrunc=diffChoice(~isnan(diffChoice));
diffChoiceTrunc=diffChoiceTrunc(1:length(diffTone(~isnan(diffTone))));
diffFullTrunc=diffFull(~isnan(diffFull));
diffFullTrunc=diffFullTrunc(1:length(diffTone(~isnan(diffTone))));
catICDiff=[diffFullTrunc diffTone(~isnan(diffTone)) diffChoiceTrunc];
subplot(1,2,2)
qqq=bar([nanmean(diffFullTrunc) nanmean(diffTone) nanmean(diffChoiceTrunc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffFull,1),1), ...
    diffFull,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffTone,1),1), ...
    diffTone,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffChoice,1),1), ...
    diffChoice,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title(['IC Inactivation']);
xticklabels({'Full Trial','Tone','Choice'});
xlabel('Condition');

wwFig.Position(3:4)=[525 275];
ylabel('Light on - Light off');
saveas(gcf,['BySess_Difference_T_MGB_IC_PercentCorrect_Opto']);
saveas(gcf,['BySess_Difference_T_MGB_IC_PercentCorrect_Opto.png']);
aovdiffMGB=anova1(catMGBDiff)
aovdiffIC=anova1(catICDiff)

end