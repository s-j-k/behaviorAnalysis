function byAnimalHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,tempTestsOnly,reinfcolor,optocolor)

ppFig=figure(15);

xVector = [1 2 1 2 1 2 1 2 1 2 1 2];
xoVector = xVector+2;
% rpc=NaN;opc=NaN;
% for jj=2:size(allDataTestsOnly,1)
%     rpc=cat(1,rpc,allDataTestsOnly{jj,27});
%     opc=cat(1,opc,allDataTestsOnly{jj,28});
% end

clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    rhittemp=mgbTempTestsOnly{jj,10}; % full MGB
    rfatemp=mgbTempTestsOnly{jj,11};
    ohittemp=mgbTempTestsOnly{jj,12};
    ofatemp=mgbTempTestsOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,1)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4) rhit(5) rfa(5) rhit(6) rfa(6)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4) ohit(5) ofa(5) ohit(6) ofa(6)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit','fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    rhittemp=mgbTempTestsOnly{jj,10}; % tone MGB
    rfatemp=mgbTempTestsOnly{jj,11};
    ohittemp = mgbTempTestsOnly{jj,14};
    ofatemp = mgbTempTestsOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,2)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4) rhit(5) rfa(5) rhit(6) rfa(6)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4) ohit(5) ofa(5) ohit(6) ofa(6)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit','fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    rhittemp=mgbTempTestsOnly{jj,10}; % choice MGB
    rfatemp=mgbTempTestsOnly{jj,11};
    ohittemp = mgbTempTestsOnly{jj,16};
    ofatemp = mgbTempTestsOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,3)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4) rhit(5) rfa(5) rhit(6) rfa(6)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4) ohit(5) ofa(5) ohit(6) ofa(6)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit','fa','hit','fa'});

% now do the IC
xVector = [1 2 1 2 1 2 1 2];
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(tempTestsOnly,1)
    rhittemp=tempTestsOnly{jj,10}; % full IC
    rfatemp=tempTestsOnly{jj,11};
    ohittemp=tempTestsOnly{jj,12};
    ofatemp=tempTestsOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,4)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector(1:8),ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal IC Full Trial Inactivation']);
xticklabels({'hit','fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(tempTestsOnly,1)
    rhittemp=tempTestsOnly{jj,10}; % tone IC
    rfatemp=tempTestsOnly{jj,11};
    ohittemp = tempTestsOnly{jj,14};
    ofatemp = tempTestsOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,5)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector(1:8),ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['IC Tone Inactivation']);
xticklabels({'hit', 'fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(tempTestsOnly,1)
    rhittemp=tempTestsOnly{jj,10}; % choice IC
    rfatemp=tempTestsOnly{jj,11};
    ohittemp = tempTestsOnly{jj,16};
    ofatemp = tempTestsOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,6)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector(1:8),ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title([' IC Choice Inactivation']);
xticklabels({'hit','fa','hit','fa'});

ppFig.Position(3:4)=[725 475];
saveas(gcf,['By Animal T_MGB_IC_HitFARate_Opto']);
saveas(gcf,['By Animal T_MGB_IC_HitFARate_Opto.png']);


xVector = [1 2 1 2 1 2 1 2 1 2 1 2];
% CONTROL by animal hit and false alarm rate
ppFig=figure(16);

% for jj=2:size(allDataCtlOnly,1)
%     if size(allDataCtlOnly{jj,27},2)>1
%         allDataCtlOnly{jj,27}=allDataCtlOnly{jj,27}';
%         rpc=cat(1,rpc,allDataCtlOnly{jj,27});
%     else
%         rpc=cat(1,rpc,allDataCtlOnly{jj,27});
%     end
%     if size(allDataCtlOnly{jj,28},2)>1
%         allDataCtlOnly{jj,28}=allDataCtlOnly{jj,28}';
%         opc=cat(1,opc,allDataCtlOnly{jj,28});
%     else
%         opc=cat(1,opc,allDataCtlOnly{jj,28});
%     end
% end

clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; % full MGB
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,12};
    ofatemp=allDataCtlOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,1)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4) rhit(5) rfa(5) rhit(6) rfa(6)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4) ohit(5) ofa(5) ohit(6) ofa(6)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit', 'fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; % tone MGB
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp = allDataCtlOnly{jj,14};
    ofatemp = allDataCtlOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,2)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4) rhit(5) rfa(5) rhit(6) rfa(6)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4) ohit(5) ofa(5) ohit(6) ofa(6)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit','fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; % choice MGB
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp = allDataCtlOnly{jj,16};
    ofatemp = allDataCtlOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,3)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4) rhit(5) rfa(5) rhit(6) rfa(6)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4) ohit(5) ofa(5) ohit(6) ofa(6)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'fa','hit','fa'});


xVector = [1 2 1 2 1 2 1 2];
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,18}; % full IC
    rfatemp=allDataCtlOnly{jj,19};
    ohittemp=allDataCtlOnly{jj,20};
    ofatemp=allDataCtlOnly{jj,21};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,4)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector(1:8),ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal IC Full Trial Inactivation']);
xticklabels({'hit', 'fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,18}; % tone IC
    rfatemp=allDataCtlOnly{jj,19};
    ohittemp = allDataCtlOnly{jj,22};
    ofatemp = allDataCtlOnly{jj,23};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,5)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector(1:8),ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['IC Tone Inactivation']);
xticklabels({'hit', 'fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,18}; % choice IC
    rfatemp=allDataCtlOnly{jj,19};
    ohittemp = allDataCtlOnly{jj,24};
    ofatemp = allDataCtlOnly{jj,25};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(2,3,6)
ppp=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(rfa,2),2), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(3),size(ohit,2),1), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3) ...
    rhit(4) rfa(4) ];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3) ...
    ohit(4) ofa(4) ];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector(1:8),ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title([' IC Choice Inactivation']);
xticklabels({'hit', 'fa','hit','fa'});

ppFig.Position(3:4)=[725 475];
saveas(gcf,['By Animal C_MGB_HitFARate_Opto']);
saveas(gcf,['By Animal C_MGB_HitFARate_Opto.png']);

end