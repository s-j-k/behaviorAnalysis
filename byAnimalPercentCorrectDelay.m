function [allDataTestsOnly]=byAnimalPercentCorrectDelay(allDataTestsOnly,firstIdxAnimal,reinfcolor,optocolor)

% test Delay 1
ee=1;
for ee=1:length(firstIdxAnimal)
    rpcVals=allDataTestsOnly{firstIdxAnimal(ee),29};
    opcVals=allDataTestsOnly{firstIdxAnimal(ee),30};
    rpc(ee)=nanmean(rpcVals);
    opc(ee)=nanmean(opcVals);
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
ee=1;
for ee=1:length(firstIdxAnimal) % Delay 2
    rpcVals=allDataTestsOnly{firstIdxAnimal(ee),31};
    opcVals=allDataTestsOnly{firstIdxAnimal(ee),32};
    rpc(ee)=nanmean(rpcVals);
    opc(ee)=nanmean(opcVals);
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
ee=1;
for ee=1:length(firstIdxAnimal) %Delay 3
    rpcVals=allDataTestsOnly{firstIdxAnimal(ee),33};
    opcVals=allDataTestsOnly{firstIdxAnimal(ee),34};
    rpc(ee)=nanmean(rpcVals);
    opc(ee)=nanmean(opcVals);
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
ee=1;
for ee=1:length(firstIdxAnimal) %Delay 4
    rpcVals=allDataTestsOnly{firstIdxAnimal(ee),35};
    opcVals=allDataTestsOnly{firstIdxAnimal(ee),36};
    rpc(ee)=nanmean(rpcVals);
    opc(ee)=nanmean(opcVals);
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
ee=1;
for ee=1:length(firstIdxAnimal) % Delay 5
    rpcVals=allDataTestsOnly{firstIdxAnimal(ee),37};
    opcVals=allDataTestsOnly{firstIdxAnimal(ee),38};
    rpc(ee)=nanmean(rpcVals);
    opc(ee)=nanmean(opcVals);
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

wwFig.Position(3:4)=[725 300];
saveas(gcf,['ByAnimal_T_MGB_Dead_PercentCorrect_Opto']);
saveas(gcf,['ByAnimal_T_MGB_Dead_PercentCorrect_Opto.png']);
end