function [mgbTempTestsOnly,allDataTestsOnly,allDataCtlOnly]=byAnimalPercentCorrect(allDataTestsOnly,allDataCtlOnly,icDataTestsOnly,reinfcolor,optocolor)

% reinfcolor1= [0.2,0.2,0.2]; %this is to plot each animal as a different
% color
% reinfcolor2= [0.45,0.45,0.45];
% reinfcolor3= [0.7,0.7,0.7];
% optocolor1=[0/255 89/255 178/255];
% optocolor2=[65/255 161/255 255/255];
% optocolor3=[181/255 216/255 255/255];

% mgbTempTestsOnly=allDataTestsOnly(1:4,:);
% mgbTempTestsOnly(5:7,:)=allDataTestsOnly(6:8,:);
mgbTempTestsOnly=allDataTestsOnly;
xVector = [1 2 1 2 1 2 1 2 1 2 1 2];
% discount sk194, fiber missing from right side
% for 194-196, the height of the implant over the MGB (which connects to
% the patch cord) was too short, so the transmittance was poor
for jj=2:size(mgbTempTestsOnly,1) % test MGB Full
    rpc(jj-1)=nanmean(mgbTempTestsOnly{jj,27});
    opc(jj-1)=nanmean(mgbTempTestsOnly{jj,28});
end
wwFig=figure(10);
subplot(2,3,1)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;

clear rpc opc 
for jj=2:size(mgbTempTestsOnly,1) % MGB Tone
    rpc(jj-1)=nanmean(mgbTempTestsOnly{jj,29});
    opc(jj-1)=nanmean(mgbTempTestsOnly{jj,30});
end
subplot(2,3,2)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc 
for jj=2:size(mgbTempTestsOnly,1) % MGB Choice
    rpc(jj-1)=nanmean(mgbTempTestsOnly{jj,31});
    opc(jj-1)=nanmean(mgbTempTestsOnly{jj,32});
end
subplot(2,3,3)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});
wwFig.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_T_MGB_PercentCorrect_Opto']);
saveas(gcf,['ByAnimal_T_MGB_PercentCorrect_Opto.png']);

% 
% wwFig=figure(15);
% clear rpc opc
% xVector = [1 2 1 2 1 2 1 2];
% for jj=2:size(icDataTestsOnly,1) % test IC Full
%     rpc(jj-1)=nanmean(icDataTestsOnly{jj,33});
%     opc(jj-1)=nanmean(icDataTestsOnly{jj,34});
% end
% % for jj=2:size(allDataTestsOnly,1) % test IC Full
% % % use only sk175 and sk195
% % tempTestsOnly(1,:)=allDataTestsOnly(1,:);
% % tempTestsOnly(2,:)=allDataTestsOnly(3,:);
% % tempTestsOnly(3:5,:)=allDataTestsOnly(5:7,:);
% % for jj=2:size(tempTestsOnly,1)
% %     rpc(jj-1)=nanmean(tempTestsOnly{jj,33});
% %     opc(jj-1)=nanmean(tempTestsOnly{jj,34});
% % end
% subplot(2,3,4)
% qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
% percentCorrect = [rpc(1) opc(1)];
% qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
% scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
%     rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
% scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
%     opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
% [h,p,ci,stats] = ttest2(rpc,opc);
% sigstar({[1,2]}, p)
% ylabel('percent correct');
% title(['IC Full Trial Inactivation']);
% xticklabels({'light off', 'light on'});
% allDataTestsOnly{jj,27}=rpc;
% allDataTestsOnly{jj,28}=opc;
% 
% clear rpc opc 
% for jj=2:size(icDataTestsOnly,1) % IC Tone
%     rpc(jj-1)=nanmean(icDataTestsOnly{jj,35});
%     opc(jj-1)=nanmean(icDataTestsOnly{jj,36});
% end
% subplot(2,3,5)
% qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
% percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
%     rpc(4) opc(4)];
% qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
% scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
%     rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
% scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
%     opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
% [h,p,ci,stats] = ttest2(rpc,opc);
% sigstar({[1,2]}, p)
% title(['IC Tone Inactivation']);
% xticklabels({'light off', 'light on'});
% 
% clear rpc opc 
% for jj=2:size(tempTestsOnly,1) % IC Choice
%     rpc(jj-1)=nanmean(tempTestsOnly{jj,37});
%     opc(jj-1)=nanmean(tempTestsOnly{jj,38});
% end
% subplot(2,3,6)
% qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
% percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
%     rpc(4) opc(4)];
% qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
% scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
%     rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
% scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
%     opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
% line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
% [h,p,ci,stats] = ttest2(rpc,opc);
% sigstar({[1,2]}, p)
% title(['IC Choice Inactivation']);
% xticklabels({'light off', 'light on'});
% 
% wwFig.Position(3:4)=[725 475];
% saveas(gcf,['ByAnimal_T_IC_PercentCorrect_Opto']);
% saveas(gcf,['ByAnimal_T_IC_PercentCorrect_Opto.png']);


clear rpc opc % Controls
xVector = [1 2 1 2 1 2 1 2 1 2 1 2];
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,27});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,28});
end
wwFig=figure(11);
subplot(2,3,1)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
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
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
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
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});


xVector = [1 2 1 2 1 2 1 2];
clear rpc opc % IC Full Trial
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,33});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,34});
end
subplot(2,3,4)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
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
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
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
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([xVector],percentCorrect, 'LineWidth', 0.5, 'Color', [0 0 0]);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['IC Choice Inactivation']);
xticklabels({'light off', 'light on'});

wwFig.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_C_MGB_IC_PercentCorrect_Opto']);
saveas(gcf,['ByAnimal_C_MGB_IC_PercentCorrect_Opto.png']);
end