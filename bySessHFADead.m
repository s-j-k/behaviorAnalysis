function bySessHFADead(allDataTestsOnly,reinfcolor,optocolor)
% By session H and FA TEST
ppFig=figure(12); % for MGB
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(1,:); % dead 1
    rfatemp=allDataTestsOnly{jj,18}(1,:);
    ohittemp=allDataTestsOnly{jj,4};
    ofatemp=allDataTestsOnly{jj,5};
    rhittemp(isnan(ohittemp))=nan;
    rfatemp(isnan(ofatemp))=nan;
    rhit=cat(1,rhit,rhittemp');
    rfa=cat(1,rfa,rfatemp');
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
title(['By Session MGB Dead 1']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(2,:); % dead 2
    rfatemp=allDataTestsOnly{jj,18}(2,:);
    ohittemp=allDataTestsOnly{jj,6};
    ofatemp=allDataTestsOnly{jj,7};
    rhittemp(isnan(ohittemp))=nan;
    rfatemp(isnan(ofatemp))=nan;
    rhit=cat(1,rhit,rhittemp');
    rfa=cat(1,rfa,rfatemp');
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
title(['MGB Dead 2']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(3,:); % dead 3
    rfatemp=allDataTestsOnly{jj,18}(3,:);
    ohittemp=allDataTestsOnly{jj,8};
    ofatemp=allDataTestsOnly{jj,9};
    rhittemp(isnan(ohittemp))=nan; % isnan should be indexing -temp vars
    rfatemp(isnan(ofatemp))=nan;
    rhit=cat(1,rhit,rhittemp');
    rfa=cat(1,rfa,rfatemp');
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
title(['By Session MGB Dead 3']);
xticklabels({'hit', 'hit','fa','fa'});

%IC 
clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(4,:); % dead 4
    rfatemp=allDataTestsOnly{jj,18}(4,:);
    ohittemp=allDataTestsOnly{jj,10};
    ofatemp=allDataTestsOnly{jj,11};
    rhittemp(isnan(ohittemp))=nan;
    rfatemp(isnan(ofatemp))=nan;
    rhit=cat(1,rhit,rhittemp');
    rfa=cat(1,rfa,rfatemp');
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
title(['By Session MGB Dead 4']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,17}(5,:); % dead 5
    rfatemp=allDataTestsOnly{jj,18}(5,:);
    ohittemp=allDataTestsOnly{jj,12};
    ofatemp=allDataTestsOnly{jj,13};
    rhittemp(isnan(ohittemp))=nan;
    rfatemp(isnan(ofatemp))=nan;
    rhit=cat(1,rhit,rhittemp');
    rfa=cat(1,rfa,rfatemp');
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
title(['MGB Dead 5']);
xticklabels({'hit', 'hit','fa','fa'});
ppFig.Position(3:4)=[725 475];
saveas(gcf,['By Sess T_MGB_Dead_HitFARate_Opto']);
saveas(gcf,['By Sess T_MGB_Dead_HitFARate_Opto.png']);
end