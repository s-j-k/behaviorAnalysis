function plotACExpertData(expertACSummaryData,reinfcolor,optocolor)

figure(105);
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
rpc=expertACSummaryData{2,2};
opc=expertACSummaryData{2,3};
rhit=expertACSummaryData{3,2}; 
rfa=expertACSummaryData{4,2}; 
ohit=expertACSummaryData{3,3}; 
ofa=expertACSummaryData{4,3}; 
rhitLat=expertACSummaryData{5,2}; 
ohitLat=expertACSummaryData{5,3}; 
rfaLat=expertACSummaryData{6,2};
ofaLat=expertACSummaryData{6,3};
rhitRate=expertACSummaryData{7,2};
ohitRate=expertACSummaryData{7,3};
rfaRate=expertACSummaryData{8,2};
ofaRate=expertACSummaryData{8,3};

% these variables are for all of the data
% need to separate out just the AC opto days, then average across animal
% then also average across session


subplot(2,2,1)
ppp=bar([nanmean(rpc) nanmean(opc)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor;];
scatter(repmat(ppp(1).XEndPoints(1),size(rpc,2),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(opc,2),2), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2],[3,4]}, [pHit])
ylabel('percent correct');
title(['By Animal AC Full Trial Inactivation']);
xticklabels({'light on','light off'});

subplot(2,2,2)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,2),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal AC Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});


subplot(2,2,3)
ppp=bar([nanmean(rhitLat) nanmean(ohitLat) nanmean(rfaLat) nanmean(ofaLat)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhitLat,2),1), ...
    rhitLat,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohitLat,2),2), ...
    ohitLat,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfaLat,2),1), ...
    rfaLat,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofaLat,2),2), ...
    ofaLat,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhitLat,ohitLat);
[h,pFA,ci,stats] = ttest2(rfaLat,ofaLat);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('Lick Latency');
title(['By Animal AC Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});


subplot(2,2,3)
ppp=bar([nanmean(rhitRate) nanmean(ohitRate) nanmean(rfaRate) nanmean(ofaRate)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhitRate,2),1), ...
    rhitRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohitRate,2),2), ...
    ohitRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfaRate,2),1), ...
    rfaRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofaRate,2),2), ...
    ofaRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhitRate,ohitRate);
[h,pFA,ci,stats] = ttest2(rfaRate,ofaRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('Lick Rate');
title(['By Animal AC Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

end
