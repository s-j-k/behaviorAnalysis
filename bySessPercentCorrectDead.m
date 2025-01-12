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

for jj=2:size(allDataTestsOnly,1) % Dead 2
    if size(allDataTestsOnly{jj,21},2)>1
        allDataTestsOnly{jj,21}=allDataTestsOnly{jj,21}';
        rpc=cat(1,rpc,allDataTestsOnly{jj,21});
    else
        rpc=cat(1,rpc,allDataTestsOnly{jj,21});
    end
    if size(allDataTestsOnly{jj,22},2)>1
        allDataTestsOnly{jj,22}=allDataTestsOnly{jj,22}';
        opc=cat(1,opc,allDataTestsOnly{jj,22});
    else
        opc=cat(1,opc,allDataTestsOnly{jj,22});
    end
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
for jj=2:size(allDataTestsOnly,1) % Dead 3
    if size(allDataTestsOnly{jj,23},2)>1
        allDataTestsOnly{jj,23}=allDataTestsOnly{jj,23}';
        rpc=cat(1,rpc,allDataTestsOnly{jj,23});
    else
        rpc=cat(1,rpc,allDataTestsOnly{jj,23});
    end
    if size(allDataTestsOnly{jj,24},2)>1
        allDataTestsOnly{jj,24}=allDataTestsOnly{jj,24}';
        opc=cat(1,opc,allDataTestsOnly{jj,24});
    else
        opc=cat(1,opc,allDataTestsOnly{jj,24});
    end
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
for jj=2:size(allDataTestsOnly,1) % Dead 4
    if size(allDataTestsOnly{jj,25},2)>1
        allDataTestsOnly{jj,25}=allDataTestsOnly{jj,25}';
        rpc=cat(1,rpc,allDataTestsOnly{jj,25});
    else
        rpc=cat(1,rpc,allDataTestsOnly{jj,25});
    end
    if size(allDataTestsOnly{jj,26},2)>1
        allDataTestsOnly{jj,26}=allDataTestsOnly{jj,26}';
        opc=cat(1,opc,allDataTestsOnly{jj,26});
    else
        opc=cat(1,opc,allDataTestsOnly{jj,26});
    end
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

clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1) % Dead 5
    if size(allDataTestsOnly{jj,27},2)>1
        allDataTestsOnly{jj,27}=allDataTestsOnly{jj,27}';
        rpc=cat(1,rpc,allDataTestsOnly{jj,27});
    else
        rpc=cat(1,rpc,allDataTestsOnly{jj,27});
    end
    if size(allDataTestsOnly{jj,28},2)>1
        allDataTestsOnly{jj,28}=allDataTestsOnly{jj,28}';
        opc=cat(1,opc,allDataTestsOnly{jj,28});
    else
        opc=cat(1,opc,allDataTestsOnly{jj,28});
    end
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