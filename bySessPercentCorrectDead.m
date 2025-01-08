function [allDataTestsOnly]=bySessPercentCorrectDead(allDataTestsOnly,reinfcolor,optocolor)

rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1) % Dead 1
    if size(allDataTestsOnly{jj,19},2)>1
        allDataTestsOnly{jj,19}=allDataTestsOnly{jj,19}';
        rpc=cat(1,rpc,allDataTestsOnly{jj,19});
    else
        rpc=cat(1,rpc,allDataTestsOnly{jj,19});
    end
    if size(allDataTestsOnly{jj,20},2)>1
        allDataTestsOnly{jj,20}=allDataTestsOnly{jj,20}';
        opc=cat(1,opc,allDataTestsOnly{jj,20});
    else
        opc=cat(1,opc,allDataTestsOnly{jj,20});
    end
end
wwFig=figure(12);
subplot(2,3,1)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Dead 1 Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,39}=rpc;
allDataTestsOnly{jj,40}=opc;

clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,21}); % need to fix dimensions here.. 
    opc=cat(1,opc,allDataTestsOnly{jj,22});
end
subplot(2,3,2)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Dead 2 Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc
rpc=NaN;opc=NaN; 
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,23});
    opc=cat(1,opc,allDataTestsOnly{jj,24});
end
subplot(2,3,3)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Dead 3 Inactivation']);
xticklabels({'light off', 'light on'});

rpc=NaN;opc=NaN; 
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,25});
    opc=cat(1,opc,allDataTestsOnly{jj,26});
end
subplot(2,3,4)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Dead 4 Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;

clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,27});
    opc=cat(1,opc,allDataTestsOnly{jj,28});
end
subplot(2,3,5)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Dead 5 Inactivation']);
xticklabels({'light off', 'light on'});

wwFig.Position(3:4)=[725 475];
saveas(gcf,['BySess_T_MGB_Dead_PercentCorrect_Opto']);
saveas(gcf,['BySess_T_MGB_Dead_PercentCorrect_Opto.png']);
end