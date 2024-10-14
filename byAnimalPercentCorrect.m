function [allDataTestsOnly,allDataCtlOnly]=byAnimalPercentCorrect(allDataTestsOnly,allDataCtlOnly)
for jj=2:size(allDataTestsOnly,1) % test MGB Full
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,27});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,28});
end
wwFig=figure(10);
subplot(2,3,1)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;

clear rpc opc 
for jj=2:size(allDataTestsOnly,1) % MGB Tone
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,29});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,30});
end
subplot(2,3,2)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc 
for jj=2:size(allDataTestsOnly,1) % MGB Choice
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,31});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,32});
end
subplot(2,3,3)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc
for jj=2:size(allDataTestsOnly,1) % test IC Full
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,33});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,34});
end
subplot(2,3,4)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['IC Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;

clear rpc opc 
for jj=2:size(allDataTestsOnly,1) % IC Tone
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,35});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,36});
end
subplot(2,3,5)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Tone Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc 
for jj=2:size(allDataTestsOnly,1) % IC Choice
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,37});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,38});
end
subplot(2,3,6)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Choice Inactivation']);
xticklabels({'light off', 'light on'});

wwFig.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_T_MGB_IC_PercentCorrect_Opto']);
saveas(gcf,['ByAnimal_T_MGB_IC_PercentCorrect_Opto.png']);


clear rpc opc % Controls
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,27});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,28});
end
wwFig=figure(10);
subplot(2,3,1)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataCtlOnly{jj,27}=rpc;
allDataCtlOnly{jj,28}=opc;

clear rpc opc
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,29});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,30});
end
subplot(2,3,2)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc
for jj=2:size(allDataCtlOnly,1) % MGB Choice 
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,31});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,32});
end
subplot(2,3,3)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc % IC Full Trial
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,33});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,34});
end
wwFig=figure(10);
subplot(2,3,4)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['IC Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataCtlOnly{jj,27}=rpc;
allDataCtlOnly{jj,28}=opc;

clear rpc opc % IC Tone
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,35});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,36});
end
subplot(2,3,5)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Tone Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc % IC Choice
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,37});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,38});
end
subplot(2,3,6)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Choice Inactivation']);
xticklabels({'light off', 'light on'});

wwFig.Position(3:4)=[725 275];
saveas(gcf,['ByAnimal_C_MGB_IC_PercentCorrect_Opto']);
saveas(gcf,['ByAnimal_C_MGB_IC_PercentCorrect_Opto.png']);