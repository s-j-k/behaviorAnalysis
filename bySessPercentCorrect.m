function [allDataTestsOnly,allDataCtlOnly]=bySessPercentCorrect(allDataTestsOnly,allDataCtlOnly)
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,27});
    opc=cat(1,opc,allDataTestsOnly{jj,28});
end
wwFig=figure(10);
subplot(2,3,1)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,29});
    opc=cat(1,opc,allDataTestsOnly{jj,30});
end
subplot(2,3,2)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,31});
    opc=cat(1,opc,allDataTestsOnly{jj,32});
end
subplot(2,3,3)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});

% now do the IC
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,33});
    opc=cat(1,opc,allDataTestsOnly{jj,34});
end
subplot(2,3,4)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;

clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,35});
    opc=cat(1,opc,allDataTestsOnly{jj,36});
end
subplot(2,3,5)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,37});
    opc=cat(1,opc,allDataTestsOnly{jj,38});
end
subplot(2,3,6)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Choice Inactivation']);
xticklabels({'light off', 'light on'});

wwFig.Position(3:4)=[725 275];
saveas(gcf,['BySess_T_MGB_PercentCorrect_Opto']);
saveas(gcf,['BySess_T_MGB_PercentCorrect_Opto.png']);

end