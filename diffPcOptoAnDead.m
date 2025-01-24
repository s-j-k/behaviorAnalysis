function diffPcOptoAnDead(allDataTestsOnly,optocolor)

diffDead1=NaN; diffDead2=NaN; diffDead3=NaN; diffDead4=NaN; diffDead5=NaN;
for jj=2:size(allDataTestsOnly,1)
    diffDead1=cat(1,diffDead1,nanmean(allDataTestsOnly{jj,49}));
    diffDead2=cat(1,diffDead2,nanmean(allDataTestsOnly{jj,50}));
    diffDead3=cat(1,diffDead3,nanmean(allDataTestsOnly{jj,51}));
    diffDead4=cat(1,diffDead4,nanmean(allDataTestsOnly{jj,52}));
    diffDead5=cat(1,diffDead5,nanmean(allDataTestsOnly{jj,53}));
end
catMGBDiff=[diffDead1(~isnan(diffDead1)) diffDead2(~isnan(diffDead2)) diffDead3(~isnan(diffDead3)) ...
    diffDead4(~isnan(diffDead4)) diffDead5(~isnan(diffDead5))];
wwFig=figure(101);
subplot(1,2,1); hold on;
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
title(['MGB Inactivation by Animal']);
xticklabels({'Dead 1','Dead 2','Dead 3','Dead 4','Dead 5'});
xlabel('Condition');

wwFig.Position(3:4)=[625 275];
ylabel('Light on - Light off');
saveas(gcf,['ByAn_Difference_T_MGB_Dead_PercentCorrect_Opto']);
saveas(gcf,['ByAn_Difference_T_MGB_Dead_PercentCorrect_Opto.png']);

aovdiffMGB=anova1(catMGBDiff)

end
