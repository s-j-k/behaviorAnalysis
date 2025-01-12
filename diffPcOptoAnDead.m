function diffPcOptoAnDead(allDataTestsOnly,optocolor)

diffFull=NaN; diffTone=NaN; diffChoice=NaN;
for jj=2:size(allDataTestsOnly,1)
    diffFull=cat(1,diffFull,nanmean(allDataTestsOnly{jj,39}));
    diffTone=cat(1,diffTone,nanmean(allDataTestsOnly{jj,40}));
    diffChoice=cat(1,diffChoice,nanmean(allDataTestsOnly{jj,41}));
end
catMGBDiff=[diffFull(~isnan(diffFull)) diffTone(~isnan(diffTone)) diffChoice(~isnan(diffChoice))];
wwFig=figure(101);
subplot(1,2,1); hold on;
qqq=bar([nanmean(diffFull) nanmean(diffTone) nanmean(diffChoice)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffFull,1),1), ...
    diffFull,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffTone,1),1), ...
    diffTone,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffChoice,1),1), ...
    diffChoice,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title(['MGB Inactivation by Animal']);
xticklabels({'Full Trial','Tone','Choice'});
xlabel('Condition');

diffFull=NaN; diffTone=NaN; diffChoice=NaN; %IC
for jj=2:size(allDataTestsOnly,1)
    diffFull=cat(1,diffFull,nanmean(allDataTestsOnly{jj,42}));
    diffTone=cat(1,diffTone,nanmean(allDataTestsOnly{jj,43}));
    diffChoice=cat(1,diffChoice,nanmean(allDataTestsOnly{jj,44}));
end
catICDiff=[diffFull(~isnan(diffFull)) diffTone(~isnan(diffTone)) diffChoice(~isnan(diffChoice))];
subplot(1,2,1)
qqq=bar([nanmean(diffFull) nanmean(diffTone) nanmean(diffChoice)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffFull,1),1), ...
    diffFull,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffTone,1),1), ...
    diffTone,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffChoice,1),1), ...
    diffChoice,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title(['MGB Inactivation by Animal']);
xticklabels({'Full Trial','Tone','Choice'});
xlabel('Condition');

wwFig.Position(3:4)=[625 275];
ylabel('Light on - Light off');
saveas(gcf,['ByAn_Difference_T_MGB_IC_PercentCorrect_Opto']);
saveas(gcf,['ByAn_Difference_T_MGB_IC_PercentCorrect_Opto.png']);

aovdiffMGB=anova1(catMGBDiff)
aovdiffIC=anova1(catICDiff)

end
