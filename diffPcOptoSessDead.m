function [allDataTestsOnly]=diffPcOptoSessDead(allDataTestsOnly,optocolor)

wwFig=figure(100);
clear jj
for jj=2:size(allDataTestsOnly,1)
    allDataTestsOnly{jj,49}=allDataTestsOnly{jj,28}-allDataTestsOnly{jj,27};
    allDataTestsOnly{jj,50}=allDataTestsOnly{jj,30}-allDataTestsOnly{jj,29};
    allDataTestsOnly{jj,51}=allDataTestsOnly{jj,32}-allDataTestsOnly{jj,31};
    allDataTestsOnly{jj,52}=allDataTestsOnly{jj,28}-allDataTestsOnly{jj,27};
    allDataTestsOnly{jj,53}=allDataTestsOnly{jj,30}-allDataTestsOnly{jj,29};
end
diffDead1=NaN; diffDead2=NaN; diffDead3=NaN;
for jj=2:size(allDataTestsOnly,1)
    diffDead1=cat(1,diffDead1,allDataTestsOnly{jj,49});
    diffDead2=cat(1,diffDead2,allDataTestsOnly{jj,50});
    diffDead3=cat(1,diffDead3,allDataTestsOnly{jj,51});
end
diffDead3Trunc=diffDead3(~isnan(diffDead3));
diffDead3Trunc=diffDead3Trunc(1:length(diffDead2(~isnan(diffDead2))));
diffDead1Trunc=diffDead1(~isnan(diffDead3));
diffDead1Trunc=diffDead1Trunc(1:length(diffDead2(~isnan(diffDead2))));
catMGBDiff=[diffDead1Trunc diffDead2(~isnan(diffDead2)) diffDead3Trunc];
hold on;
qqq=bar([nanmean(diffDead1) nanmean(diffDead2) nanmean(diffDead3)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffDead1,1),1), ...
    diffDead1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffDead2,1),1), ...
    diffDead2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffDead3,1),1), ...
    diffDead3,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title(['MGB Inactivation']);
xticklabels({'Full Trial','Tone','Choice'});
xlabel('Condition');

saveas(gcf,['BySess_Difference_T_MGB_Dead_PercentCorrect_Opto']);
saveas(gcf,['BySess_Difference_T_MGB_Dead_PercentCorrect_Opto.png']);
aovdiffMGB=anova1(catMGBDiff)
% aovdiffIC=anova1(catICDiff)

end