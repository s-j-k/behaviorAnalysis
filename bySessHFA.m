function bySessHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,tempTestsOnly,reinfcolor,optocolor)

% By session H and FA TEST
ppFig=figure(12); % for MGB
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    rhittemp=mgbTempTestsOnly{jj,10}; 
    rfatemp=mgbTempTestsOnly{jj,11};
    ohittemp=mgbTempTestsOnly{jj,12};
    ofatemp=mgbTempTestsOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,1)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    rhittemp=mgbTempTestsOnly{jj,10}; 
    rfatemp=mgbTempTestsOnly{jj,11};
    ohittemp=mgbTempTestsOnly{jj,14};
    ofatemp=mgbTempTestsOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,2)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(mgbTempTestsOnly,1)
    rhittemp=mgbTempTestsOnly{jj,10}; 
    rfatemp=mgbTempTestsOnly{jj,11};
    ohittemp=mgbTempTestsOnly{jj,16};
    ofatemp=mgbTempTestsOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,3)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

%IC 
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(tempTestsOnly,1)
    rhittemp=tempTestsOnly{jj,18}; 
    rfatemp=tempTestsOnly{jj,19};
    ohittemp=tempTestsOnly{jj,20};
    ofatemp=tempTestsOnly{jj,21};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,4)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal IC Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(tempTestsOnly,1)
    rhittemp=tempTestsOnly{jj,18}; 
    rfatemp=tempTestsOnly{jj,19};
    ohittemp=tempTestsOnly{jj,22};
    ofatemp=tempTestsOnly{jj,23};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,5)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['IC Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(tempTestsOnly,1)
    rhittemp=tempTestsOnly{jj,18}; 
    rfatemp=tempTestsOnly{jj,19};
    ohittemp=tempTestsOnly{jj,24};
    ofatemp=tempTestsOnly{jj,25};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,6)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' IC Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

ppFig.Position(3:4)=[725 475];
saveas(gcf,['By Sess T_MGB_IC_HitFARate_Opto']);
saveas(gcf,['By Sess T_MGB_IC_HitFARate_Opto.png']);

%now do CTLs
ppFig=figure(13);
clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; 
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,12};
    ofatemp=allDataCtlOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,1)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; 
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,14};
    ofatemp=allDataCtlOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,2)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; 
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,16};
    ofatemp=allDataCtlOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,3)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

%IC
clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,18}; 
    rfatemp=allDataCtlOnly{jj,19};
    ohittemp=allDataCtlOnly{jj,20};
    ofatemp=allDataCtlOnly{jj,21};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,4)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal IC Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,18}; 
    rfatemp=allDataCtlOnly{jj,19};
    ohittemp=allDataCtlOnly{jj,22};
    ofatemp=allDataCtlOnly{jj,23};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,5)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['IC Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,18}; 
    rfatemp=allDataCtlOnly{jj,19};
    ohittemp=allDataCtlOnly{jj,24};
    ofatemp=allDataCtlOnly{jj,25};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(2,3,6)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

ppFig.Position(3:4)=[725 475];
saveas(gcf,['By Sess C_MGB_IC_HitFARate_Opto']);
saveas(gcf,['By Sess C_MGB_IC_HitFARate_Opto.png']);

end