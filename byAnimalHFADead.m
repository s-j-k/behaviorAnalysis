function byAnimalHFADead(allDataTestsOnly,reinfcolor,optocolor)

ppFig=figure(15);
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
line([1 2],hitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
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
line([1 2],hitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
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
line([1 2],hitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
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
line([1 2],hitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
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
line([1 2],hitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],hitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([3 4],ohitAndFa(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,3],[2,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Dead 5']);
xticklabels({'hit', 'hit','fa','fa'});
ppFig.Position(3:4)=[725 475];
saveas(gcf,['By Animal T_MGB_Dead_HitFARate_Opto']);
saveas(gcf,['By Animal T_MGB_Dead_HitFARate_Opto.png']);
%% working on this
% test animals, plot d' with the bar plot, per animal 
clear dprimeFIg
dprimeFig=figure(21);
hold on;
clear rhittemp rfatemp rfatemp2 rhittemp2
for jj=2:size(allDataTestsOnly,1)
    rhittemp(:,jj-1)=allDataTestsOnly{jj,17}(1,:);
    rfatemp(:,jj-1)=allDataTestsOnly{jj,18}(1,:);
    ohittemp(:,jj-1)=allDataTestsOnly{jj,4};
    ofatemp(:,jj-1)=allDataTestsOnly{jj,5};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp(:,jj-1));
    rfa(jj-1)=nanmean(rfatemp(:,jj-1));
    ohit(jj-1)=nanmean(ohittemp(:,jj-1));
    ofa(jj-1) = nanmean(ofatemp(:,jj-1));
end
rhittemp(rhittemp==1)=1-1./(2*140);
rfatemp(rfatemp==1)=1-1./(2*140);
ohittemp(ohittemp==1)=1-1./(2*140);
ofatemp(ofatemp==1)=1-1./(2*140);
rhittemp(rhittemp==0)=1./(2*140);
rfatemp(rfatemp==0)=1./(2*140);
ohittemp(ohittemp==0)=1./(2*140);
ofatemp(ofatemp==0)=1./(2*140);
[dp_r,c_r]=dprime_simple(rhittemp,rfatemp);
[dp_opto,c_opto]=dprime_simple(ohittemp, ofatemp);
subplot(2, 3, 1); 
dp_re=nansem(dp_r');
dp_rm=nanmean(dp_r');
hold on;
odp_re=nansem(dp_opto');
odp_rm=nanmean(dp_opto');
dp_opto_An=nanmean(dp_opto);
dp_r_An=mean(dp_r);
ddd=bar([nanmean(dp_rm) nanmean(odp_rm)]); hold on;
ddd(1).FaceColor='flat';ddd(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(ddd(1).XEndPoints(1),size(mean(dp_r),2),1), ...
    mean(dp_r),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ddd(1).XEndPoints(2),size(nanmean(dp_opto),2),1), ...
    nanmean(dp_opto),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(nanmean(dp_r),nanmean(dp_opto));
sigstar({[1,2]}, p)
dprimeVals = [dp_r_An(1) dp_opto_An(1) dp_r_An(2) dp_opto_An(2) dp_r_An(3) dp_opto_An(3)];
dprimeVals=dprimeVals(~isnan(dprimeVals));
line([1 2],dprimeVals(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
box off
allDataTestsOnly{2,45}=p;
allDataTestsOnly{2,39} = [dp_r_An; dp_opto_An];
title('Delay 1');ylabel('d prime');ylim([-1 4]);
xticklabels({'','light off','light on'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp(:,jj-1)=allDataTestsOnly{jj,17}(2,:);
    rfatemp(:,jj-1)=allDataTestsOnly{jj,18}(2,:);
    ohittemp(:,jj-1)=allDataTestsOnly{jj,6};
    ofatemp(:,jj-1)=allDataTestsOnly{jj,7};
    rhittemp(isnan(ohittemp))=nan;
    rfatemp(isnan(ofatemp))=nan;
    rhit(jj-1)=nanmean(rhittemp(:,jj-1));
    rfa(jj-1)=nanmean(rfatemp(:,jj-1));
    ohit(jj-1)=nanmean(ohittemp(:,jj-1));
    ofa(jj-1) = nanmean(ofatemp(:,jj-1));
end
rhittemp(rhittemp==1)=1-1./(2*140);
rfatemp(rfatemp==1)=1-1./(2*140);
ohittemp(ohittemp==1)=1-1./(2*140);
ofatemp(ofatemp==1)=1-1./(2*140);
rhittemp(rhittemp==0)=1./(2*140);
rfatemp(rfatemp==0)=1./(2*140);
ohittemp(ohittemp==0)=1./(2*140);
ofatemp(ofatemp==0)=1./(2*140);
[dp_r,c_r]=dprime_simple(rhittemp,rfatemp);
[dp_opto,c_opto]=dprime_simple(ohittemp, ofatemp);
subplot(2, 3, 2); 
dp_re=nansem(dp_r');
dp_rm=nanmean(dp_r');
hold on;
odp_re=nansem(dp_opto');
odp_rm=nanmean(dp_opto');
dp_opto_An=nanmean(dp_opto);
dp_r_An=nanmean(dp_r);
ddd=bar([nanmean(dp_rm) nanmean(odp_rm)]); hold on;
ddd(1).FaceColor='flat';ddd(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(ddd(1).XEndPoints(1),size(mean(dp_r),2),1), ...
    mean(dp_r),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ddd(1).XEndPoints(2),size(nanmean(dp_opto),2),1), ...
    nanmean(dp_opto),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(nanmean(dp_r),nanmean(dp_opto));
sigstar({[1,2]}, p)
dprimeVals = [dp_r_An(1) dp_opto_An(1) dp_r_An(2) dp_opto_An(2) dp_r_An(3) dp_opto_An(3)];
dprimeVals=dprimeVals(~isnan(dprimeVals));
line([1 2],dprimeVals(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
allDataTestsOnly{2,45}=p;
allDataTestsOnly{2,39} = [dp_r_An; dp_opto_An];
box off
title('Delay 2');ylim([-1 4]);
xticklabels({'','','light off','','light on'});

for jj=2:size(allDataTestsOnly,1)
    rhittemp(:,jj-1)=allDataTestsOnly{jj,17}(3,:);
    rfatemp(:,jj-1)=allDataTestsOnly{jj,18}(3,:);
    ohittemp(:,jj-1)=allDataTestsOnly{jj,8};
    ofatemp(:,jj-1)=allDataTestsOnly{jj,9};
    rhittemp(isnan(ohittemp))=nan;
    rfatemp(isnan(ofatemp))=nan;
    rhit(jj-1)=nanmean(rhittemp(:,jj-1));
    rfa(jj-1)=nanmean(rfatemp(:,jj-1));
    ohit(jj-1)=nanmean(ohittemp(:,jj-1));
    ofa(jj-1) = nanmean(ofatemp(:,jj-1));
end
rhittemp(rhittemp==1)=1-1./(2*140);
rfatemp(rfatemp==1)=1-1./(2*140);
ohittemp(ohittemp==1)=1-1./(2*140);
ofatemp(ofatemp==1)=1-1./(2*140);
rhittemp(rhittemp==0)=1./(2*140);
rfatemp(rfatemp==0)=1./(2*140);
ohittemp(ohittemp==0)=1./(2*140);
ofatemp(ofatemp==0)=1./(2*140);
[dp_r,c_r]=dprime_simple(rhittemp,rfatemp);
[dp_opto,c_opto]=dprime_simple(ohittemp, ofatemp);
subplot(2, 3, 3); 
dp_re=nansem(dp_r');
dp_rm=nanmean(dp_r');
hold on;
odp_re=nansem(dp_opto');
odp_rm=nanmean(dp_opto');
dp_opto_An=nanmean(dp_opto);
dp_r_An=mean(dp_r);
ddd=bar([nanmean(dp_rm) nanmean(odp_rm)]); hold on;
ddd(1).FaceColor='flat';ddd(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(ddd(1).XEndPoints(1),size(mean(dp_r),2),1), ...
    mean(dp_r),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ddd(1).XEndPoints(2),size(nanmean(dp_opto),2),1), ...
    nanmean(dp_opto),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(nanmean(dp_r),nanmean(dp_opto));
sigstar({[1,2]}, p)
dprimeVals = [dp_r_An(1) dp_opto_An(1) dp_r_An(2) dp_opto_An(2) dp_r_An(3) dp_opto_An(3)];
dprimeVals=dprimeVals(~isnan(dprimeVals));
line([1 2],dprimeVals(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
allDataTestsOnly{2,45}=p;
allDataTestsOnly{2,39} = [dp_r_An; dp_opto_An];
box off
title('Delay 3');ylim([-1 4]);
xticklabels({'','light off','light on'});

for jj=2:size(allDataTestsOnly,1)
    rhittemp(:,jj-1)=allDataTestsOnly{jj,17}(4,:);
    rfatemp(:,jj-1)=allDataTestsOnly{jj,18}(4,:);
    ohittemp(:,jj-1)=allDataTestsOnly{jj,10};
    ofatemp(:,jj-1)=allDataTestsOnly{jj,11};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp(:,jj-1));
    rfa(jj-1)=nanmean(rfatemp(:,jj-1));
    ohit(jj-1)=nanmean(ohittemp(:,jj-1));
    ofa(jj-1) = nanmean(ofatemp(:,jj-1));
end
rhittemp(rhittemp==1)=1-1./(2*140);
rfatemp(rfatemp==1)=1-1./(2*140);
ohittemp(ohittemp==1)=1-1./(2*140);
ofatemp(ofatemp==1)=1-1./(2*140);
rhittemp(rhittemp==0)=1./(2*140);
rfatemp(rfatemp==0)=1./(2*140);
ohittemp(ohittemp==0)=1./(2*140);
ofatemp(ofatemp==0)=1./(2*140);
[dp_r,c_r]=dprime_simple(rhittemp,rfatemp);
[dp_opto,c_opto]=dprime_simple(ohittemp, ofatemp);
subplot(2, 3, 4); 
dp_re=nansem(dp_r');
dp_rm=nanmean(dp_r');
hold on;
odp_re=nansem(dp_opto');
odp_rm=nanmean(dp_opto');
dp_opto_An=nanmean(dp_opto);
dp_r_An=mean(dp_r);
ddd=bar([nanmean(dp_rm) nanmean(odp_rm)]); hold on;
ddd(1).FaceColor='flat';ddd(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(ddd(1).XEndPoints(1),size(mean(dp_r),2),1), ...
    mean(dp_r),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ddd(1).XEndPoints(2),size(nanmean(dp_opto),2),1), ...
    nanmean(dp_opto),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(nanmean(dp_r),nanmean(dp_opto));
sigstar({[1,2]}, p)
dprimeVals = [dp_r_An(1) dp_opto_An(1) dp_r_An(2) dp_opto_An(2) dp_r_An(3) dp_opto_An(3)];
dprimeVals=dprimeVals(~isnan(dprimeVals));
line([1 2],dprimeVals(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
allDataTestsOnly{2,45}=p;
allDataTestsOnly{2,39} = [dp_r_An; dp_opto_An];
box off
title('Delay 4');ylim([-1 4]);
xticklabels({'','light off','light on'});

for jj=2:size(allDataTestsOnly,1)
    rhittemp(:,jj-1)=allDataTestsOnly{jj,17}(5,:);
    rfatemp(:,jj-1)=allDataTestsOnly{jj,18}(5,:);
    ohittemp(:,jj-1)=allDataTestsOnly{jj,12};
    ofatemp(:,jj-1)=allDataTestsOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp(:,jj-1));
    rfa(jj-1)=nanmean(rfatemp(:,jj-1));
    ohit(jj-1)=nanmean(ohittemp(:,jj-1));
    ofa(jj-1) = nanmean(ofatemp(:,jj-1));
end
rhittemp(rhittemp==1)=1-1./(2*140);
rfatemp(rfatemp==1)=1-1./(2*140);
ohittemp(ohittemp==1)=1-1./(2*140);
ofatemp(ofatemp==1)=1-1./(2*140);
rhittemp(rhittemp==0)=1./(2*140);
rfatemp(rfatemp==0)=1./(2*140);
ohittemp(ohittemp==0)=1./(2*140);
ofatemp(ofatemp==0)=1./(2*140);
[dp_r,c_r]=dprime_simple(rhittemp,rfatemp);
[dp_opto,c_opto]=dprime_simple(ohittemp, ofatemp);
subplot(2,3,5); 
dp_re=nansem(dp_r');
dp_rm=nanmean(dp_r');
hold on;
odp_re=nansem(dp_opto');
odp_rm=nanmean(dp_opto');
dp_opto_An=nanmean(dp_opto);
dp_r_An=mean(dp_r);
ddd=bar([nanmean(dp_rm) nanmean(odp_rm)]); hold on;
ddd(1).FaceColor='flat';ddd(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(ddd(1).XEndPoints(1),size(mean(dp_r),2),1), ...
    mean(dp_r),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ddd(1).XEndPoints(2),size(nanmean(dp_opto),2),1), ...
    nanmean(dp_opto),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(nanmean(dp_r),nanmean(dp_opto));
sigstar({[1,2]}, p)
dprimeVals = [dp_r_An(1) dp_opto_An(1) dp_r_An(2) dp_opto_An(2) dp_r_An(3) dp_opto_An(3)];
dprimeVals=dprimeVals(~isnan(dprimeVals));
line([1 2],dprimeVals(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
line([1 2],dprimeVals(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
allDataTestsOnly{2,45}=p;
allDataTestsOnly{2,39} = [dp_r_An; dp_opto_An];
box off
title('Delay 5');ylim([-1 4]);
xticklabels({'','light off','light on'});


dprimeFig.Position(3:4)=[725 525];
saveas(gcf,['By Animal_T_MGB_dp_bar_Delay_Opto']);
saveas(gcf,['By Animal_T_MGB_dp_bar_Delay_Opto.png']);
saveas(gcf,['By Animal_T_MGB_dp_bar_Delay_Opto.pdf']);

end