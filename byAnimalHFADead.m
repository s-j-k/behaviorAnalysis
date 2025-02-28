function byAnimalHFADead(allDataTestsOnly,reinfcolor,optocolor)

ppFig=figure(15);
xVector = [1 2 1 2 1 2];
xoVector = xVector+2;
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(1,:); % dead 1
    rfatemp=allDataTestsOnly{jj,18}(1,:);
    ohittemp=allDataTestsOnly{jj,4};
    ofatemp=allDataTestsOnly{jj,5};
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
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);

[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Dead 1']);
xticklabels({'hit', 'fa','hit','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(2,:); % dead 2
    rfatemp=allDataTestsOnly{jj,18}(2,:);
    ohittemp=allDataTestsOnly{jj,6};
    ofatemp=allDataTestsOnly{jj,7};
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
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);

[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Dead 2']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(3,:); % dead 3
    rfatemp=allDataTestsOnly{jj,18}(3,:);
    ohittemp=allDataTestsOnly{jj,8};
    ofatemp=allDataTestsOnly{jj,9};
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
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);

[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Dead 3']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(4,:); % dead 4
    rfatemp=allDataTestsOnly{jj,18}(4,:);
    ohittemp=allDataTestsOnly{jj,10};
    ofatemp=allDataTestsOnly{jj,11};
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
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);

[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Dead 4']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(5,:); % dead 5
    rfatemp=allDataTestsOnly{jj,18}(5,:);
    ohittemp=allDataTestsOnly{jj,12};
    ofatemp=allDataTestsOnly{jj,13};
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
hitAndFa = [rhit(1) rfa(1) rhit(2) rfa(2) rhit(3) rfa(3)];
ohitAndFa = [ohit(1) ofa(1) ohit(2) ofa(2) ohit(3) ofa(3)];
line(xVector,hitAndFa, 'LineWidth', 0.5, 'Color', [0 0 0]);
line(xoVector,ohitAndFa,'LineWidth', 0.5, 'Color', [0 0 0]);

[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Dead 5']);
xticklabels({'hit', 'hit','fa','fa'});
ppFig.Position(3:4)=[725 475];
saveas(gcf,['By Animal T_MGB_Dead_HitFARate_Opto']);
saveas(gcf,['By Animal T_MGB_Dead_HitFARate_Opto.png']);
end