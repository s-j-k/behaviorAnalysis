function [allDataTestsOnly]=byAnimalPercentCorrectDead(allDataTestsOnly,reinfcolor,optocolor)

for jj=2:size(allDataTestsOnly,1) % test Dead 1
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,19});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,20});
    allDataTestsOnly{jj,29}=rpc(jj-1);
    allDataTestsOnly{jj,30}=opc(jj-1);
end
wwFig=figure(10);
subplot(2,3,1)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Dead 1 Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc 
for jj=2:size(allDataTestsOnly,1) % MGB Dead 2
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,21});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,22});
    allDataTestsOnly{jj,31}=rpc(jj-1);
    allDataTestsOnly{jj,32}=opc(jj-1);
end
subplot(2,3,2)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Dead 2 Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc 
for jj=2:size(allDataTestsOnly,1) % MGB Dead 3
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,23});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,24});
    allDataTestsOnly{jj,33}=rpc(jj-1);
    allDataTestsOnly{jj,34}=opc(jj-1);
end
subplot(2,3,3)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Dead 3 Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc
for jj=2:size(allDataTestsOnly,1) % MGB Dead 4
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,25});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,26});
    allDataTestsOnly{jj,35}=rpc(jj-1);
    allDataTestsOnly{jj,36}=opc(jj-1);
end
subplot(2,3,4)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Dead 4 Inactivation']);
xticklabels({'light off', 'light on'});

clear rpc opc 
for jj=2:size(allDataTestsOnly,1) % MGB Dead 5
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,27});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,28});
    allDataTestsOnly{jj,37}=rpc(jj-1);
    allDataTestsOnly{jj,38}=opc(jj-1);
end
subplot(2,3,5)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3)];
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Dead 5 Inactivation']);
xticklabels({'light off', 'light on'});

wwFig.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_T_MGB_Dead_PercentCorrect_Opto']);
saveas(gcf,['ByAnimal_T_MGB_Dead_PercentCorrect_Opto.png']);
end