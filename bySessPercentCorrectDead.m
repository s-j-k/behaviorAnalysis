function [allDataTestsOnly,allDataCtlOnly]=bySessPercentCorrectDead(allDataTestsOnly,allDataCtlOnly,tempTestsOnly,mgbTempTestsOnly,reinfcolor,optocolor)

% reinfcolor1= [0.2,0.2,0.2]; %this is to plot each animal as a different
% color
% reinfcolor2= [0.45,0.45,0.45];
% reinfcolor3= [0.7,0.7,0.7];
% optocolor1=[0/255 89/255 178/255];
% optocolor2=[65/255 161/255 255/255];
% optocolor3=[181/255 216/255 255/255];

rpc=NaN;opc=NaN;
for jj=2:size(mgbTempTestsOnly,1) % test MGB Full
    if size(mgbTempTestsOnly{jj,27},2)>1
        mgbTempTestsOnly{jj,27}=mgbTempTestsOnly{jj,27}';
        rpc=cat(1,rpc,mgbTempTestsOnly{jj,27});
    else
        rpc=cat(1,rpc,mgbTempTestsOnly{jj,27});
    end
    if size(mgbTempTestsOnly{jj,28},2)>1
        mgbTempTestsOnly{jj,28}=mgbTempTestsOnly{jj,28}';
        opc=cat(1,opc,mgbTempTestsOnly{jj,28});
    else
        opc=cat(1,opc,mgbTempTestsOnly{jj,28});
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
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    rpc=cat(1,rpc,mgbTempTestsOnly{jj,29});
    opc=cat(1,opc,mgbTempTestsOnly{jj,30});
end
subplot(2,3,2)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN; 
for jj=2:size(mgbTempTestsOnly,1)
    rpc=cat(1,rpc,mgbTempTestsOnly{jj,31});
    opc=cat(1,opc,mgbTempTestsOnly{jj,32});
end
subplot(2,3,3)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});

% now do the IC
rpc=NaN;opc=NaN; 
for jj=2:size(tempTestsOnly,1)% IC Full
    rpc=cat(1,rpc,tempTestsOnly{jj,33});
    opc=cat(1,opc,tempTestsOnly{jj,34});
end
subplot(2,3,4)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['IC Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;

clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(tempTestsOnly,1) % IC Tone
    rpc=cat(1,rpc,tempTestsOnly{jj,35});
    opc=cat(1,opc,tempTestsOnly{jj,36});
end
subplot(2,3,5)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(tempTestsOnly,1) % IC Choice
    rpc=cat(1,rpc,tempTestsOnly{jj,37});
    opc=cat(1,opc,tempTestsOnly{jj,38});
end
subplot(2,3,6)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
% scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
%     rpc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
% scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
%     opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Choice Inactivation']);
xticklabels({'light off', 'light on'});

wwFig.Position(3:4)=[725 475];
saveas(gcf,['BySess_T_MGB_IC_PercentCorrect_Opto']);
saveas(gcf,['BySess_T_MGB_IC_PercentCorrect_Opto.png']);


% now do the controls
allDataCtlOnly{5,27}=allDataCtlOnly{5,27}';
allDataCtlOnly{5,28}=allDataCtlOnly{5,28}';
clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    if size(allDataCtlOnly{jj,27},2)>1
        allDataCtlOnly{jj,27}=allDataCtlOnly{jj,27}';
        rpc=cat(1,rpc,allDataCtlOnly{jj,27});
    else
        rpc=cat(1,rpc,allDataCtlOnly{jj,27});
    end
    if size(allDataCtlOnly{jj,28},2)>1
        allDataCtlOnly{jj,28}=allDataCtlOnly{jj,28}';
        opc=cat(1,opc,allDataCtlOnly{jj,28});
    else
        opc=cat(1,opc,allDataCtlOnly{jj,28});
    end
end
wwFig=figure(13);
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
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataCtlOnly{jj,27}=rpc;
allDataCtlOnly{jj,28}=opc;


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,29});
    opc=cat(1,opc,allDataCtlOnly{jj,30});
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
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,31});
    opc=cat(1,opc,allDataCtlOnly{jj,32});
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
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,33});
    opc=cat(1,opc,allDataCtlOnly{jj,34});
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
title(['IC Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataCtlOnly{jj,27}=rpc;
allDataCtlOnly{jj,28}=opc;


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,35});
    opc=cat(1,opc,allDataCtlOnly{jj,36});
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
title(['IC Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,37});
    opc=cat(1,opc,allDataCtlOnly{jj,38});
end
subplot(2,3,6)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
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
saveas(gcf,['BySess_C_MGB_IC_PercentCorrect_Opto']);
saveas(gcf,['BySess_C_MGB_IC_PercentCorrect_Opto.png']);

end