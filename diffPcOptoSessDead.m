function [allDataTestsOnly]=diffPcOptoSessDead(allDataTestsOnly,optocolor)
wwFig=figure(100);
clear jj
for jj=2:size(allDataTestsOnly,1)
    allDataTestsOnly{jj,49}=allDataTestsOnly{jj,40}-allDataTestsOnly{jj,39};
    allDataTestsOnly{jj,50}=allDataTestsOnly{jj,42}-allDataTestsOnly{jj,41};
    allDataTestsOnly{jj,51}=allDataTestsOnly{jj,44}-allDataTestsOnly{jj,43};
    allDataTestsOnly{jj,52}=allDataTestsOnly{jj,46}-allDataTestsOnly{jj,45};
    allDataTestsOnly{jj,53}=allDataTestsOnly{jj,48}-allDataTestsOnly{jj,47};
end
diffDead1=NaN; diffDead2=NaN; diffDead3=NaN; diffDead4=NaN; diffDead5 = NaN;
for jj=2:size(allDataTestsOnly,1)
    diffDead1=cat(1,diffDead1,allDataTestsOnly{jj,49});
    diffDead2=cat(1,diffDead2,allDataTestsOnly{jj,50});
    diffDead3=cat(1,diffDead3,allDataTestsOnly{jj,51});
    diffDead4 =cat(1,diffDead4,allDataTestsOnly{jj,52});
    diffDead5 =cat(1,diffDead5,allDataTestsOnly{jj,53});
end
catMGBDiff=[diffDead1 diffDead2 diffDead3 diffDead4 diffDead5];
hold on;
qqq=bar([nanmean(diffDead1) nanmean(diffDead2) nanmean(diffDead3) ...
    nanmean(diffDead4) nanmean(diffDead5)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffDead1,1),1), ...
    diffDead1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffDead2,1),1), ...
    diffDead2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffDead3,1),1), ...
    diffDead3,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(4),size(diffDead4,1),1), ...
    diffDead4,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(5),size(diffDead5,1),1), ...
    diffDead5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title(['MGB Inactivation']);
xticklabels({'','Dead 1','Dead 2','Dead 3', 'Dead 4','Dead 5'});
xlabel('Condition');

saveas(gcf,['BySess_Difference_T_MGB_Dead_PercentCorrect_Opto']);
saveas(gcf,['BySess_Difference_T_MGB_Dead_PercentCorrect_Opto.png']);
aovdiffMGB=anova1(catMGBDiff)
% aovdiffIC=anova1(catICDiff)

end