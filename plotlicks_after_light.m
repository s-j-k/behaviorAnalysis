function plotlicks_after_light(lickMat,allDataTestsOnly, allDataCtlOnly,icDataTestsOnly)

% for sk198 39+42 - full and choice, 40+44 - tone and choice

nbsubj=1;
lickhistr_miss = lickhistrmiss{nbsubj};
lickhistr_cr = lickhistrcr{nbsubj};
lickhisto_miss= lickhistomiss{nbsubj};
lickhisto_cr = lickhistocr{nbsubj};
lickhistot_miss= lickhistotmiss{nbsubj};
lickhistot_cr=lickhistotcr{nbsubj};
lickhistoc_miss=lickhistocmiss{nbsubj};
lickhistoc_cr=lickhistoccr{nbsubj};

figlight=figure; title('sk198');
subplot(2,4,1);title('full trial, miss, d39');hold on;
plot(bins,lickhisto_miss(39,:),'b');
plot(bins,lickhistr_miss(39,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,2);title('full trial, miss, d42');hold on;
plot(bins,lickhisto_miss(42,:),'b');
plot(bins,lickhistr_miss(42,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,3);title('choice, miss, d39');hold on;
plot(bins,lickhistoc_miss(39,:),'b');
plot(bins,lickhistr_miss(39,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,4);title('choice, miss, d40');hold on;
plot(bins,lickhistoc_miss(40,:),'b');
plot(bins,lickhistr_miss(40,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,5);title('tone, miss, d40');hold on;
plot(bins,lickhistot_miss(40,:),'b');
plot(bins,lickhistr_miss(40,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,6);title('tone, miss, d44');hold on;
plot(bins,lickhistot_miss(44,:),'b');
plot(bins,lickhistr_miss(44,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,7);title('choice, miss, d42');hold on;
plot(bins,lickhistoc_miss(42,:),'b');
plot(bins,lickhistr_miss(42,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,8);title('choice, miss, d44');hold on;
plot(bins,lickhistoc_miss(44,:),'b');
plot(bins,lickhistr_miss(44,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);
filetype='fig';
saveas(figlight,[subjlist{nbsubj} '_licks_after_light' num2str(nbins) 'bins.' filetype]);

%% sk198
lickhistr_hit =lickhistrhit{nbsubj};
lickhistr_fa =lickhistrfa{nbsubj};
lickhisto_hit =lickhistohit{nbsubj};
lickhisto_fa= lickhistofa{nbsubj};
lickhistot_hit = lickhistothit{nbsubj};
lickhistot_fa = lickhistotfa{nbsubj};
lickhistoc_hit = lickhistochit{nbsubj};
lickhistoc_fa = lickhistocfa{nbsubj};

figAction=figure; title('sk198');
subplot(2,4,1);title('full trial hit d39');hold on;
plot(bins,lickhisto_hit(39,:),'b');
plot(bins,lickhistr_hit(39,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,2);title('full trial hit');hold on;
plot(bins,lickhisto_hit(42,:),'b');
plot(bins,lickhistr_hit(42,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,3);title('choice');hold on;
plot(bins,lickhistoc_hit(39,:),'b');
plot(bins,lickhistr_hit(39,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,4);title('choice');hold on;
plot(bins,lickhistoc_hit(42,:),'b');
plot(bins,lickhistr_hit(42,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,5);title('tone');hold on;
plot(bins,lickhistot_hit(40,:),'b');
plot(bins,lickhistr_hit(40,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,6);title('tone');hold on;
plot(bins,lickhistot_hit(44,:),'b');
plot(bins,lickhistr_hit(44,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,7);title('choice');hold on;
plot(bins,lickhistoc_hit(40,:),'b');
plot(bins,lickhistr_hit(40,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,8);title('choice');hold on;
plot(bins,lickhistoc_hit(44,:),'b');
plot(bins,lickhistr_hit(44,:),'k');
xline(2.8);xline(0.3);axis tight;
filetype='fig';
saveas(figAction,[subjlist{nbsubj} '_licks_after_light_Action' num2str(nbins) 'bins.' filetype]);


%% 9, 10, 14, 15 for sk203 
nbsubj=6;
lickhistr_miss = lickhistrmiss{nbsubj};
lickhistr_cr = lickhistrcr{nbsubj};
lickhisto_miss= lickhistomiss{nbsubj};
lickhisto_cr = lickhistocr{nbsubj};
lickhistot_miss= lickhistotmiss{nbsubj};
lickhistot_cr=lickhistotcr{nbsubj};
lickhistoc_miss=lickhistocmiss{nbsubj};
lickhistoc_cr=lickhistoccr{nbsubj};

figlight=figure;  title('sk203');
subplot(2,4,1);title('full trial');hold on;
plot(bins,lickhisto_miss(9,:),'b');
plot(bins,lickhistr_miss(9,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,2);title('full trial');hold on;
plot(bins,lickhisto_miss(14,:),'b');
plot(bins,lickhistr_miss(14,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,3);title('choice');hold on;
plot(bins,lickhistoc_miss(9,:),'b');
plot(bins,lickhistr_miss(9,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,4);title('choice');hold on;
plot(bins,lickhistoc_miss(14,:),'b');
plot(bins,lickhistr_miss(14,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,5);title('tone');hold on;
plot(bins,lickhistot_miss(10,:),'b');
plot(bins,lickhistr_miss(10,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,6);title('tone');hold on;
plot(bins,lickhistot_miss(15,:),'b');
plot(bins,lickhistr_miss(15,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,7);title('choice');hold on;
plot(bins,lickhistoc_miss(10,:),'b');
plot(bins,lickhistr_miss(10,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,8);title('choice');hold on;
plot(bins,lickhistoc_miss(15,:),'b');
plot(bins,lickhistr_miss(15,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);
filetype='fig';
saveas(figlight,[subjlist{nbsubj} '_licks_after_light' num2str(nbins) 'bins.' filetype]);

% sk203
lickhistr_hit =lickhistrhit{nbsubj};
lickhistr_fa =lickhistrfa{nbsubj};
lickhisto_hit =lickhistohit{nbsubj};
lickhisto_fa= lickhistofa{nbsubj};
lickhistot_hit = lickhistothit{nbsubj};
lickhistot_fa = lickhistotfa{nbsubj};
lickhistoc_hit = lickhistochit{nbsubj};
lickhistoc_fa = lickhistocfa{nbsubj};


figAction=figure; 
subplot(2,4,1);title('full trial hit');hold on;
plot(bins,lickhisto_hit(9,:),'b');
plot(bins,lickhistr_hit(9,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,2);title('full trial hit');hold on;
plot(bins,lickhisto_hit(14,:),'b');
plot(bins,lickhistr_hit(14,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,3);title('choice');hold on;
plot(bins,lickhistoc_hit(9,:),'b');
plot(bins,lickhistr_hit(9,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,4);title('choice');hold on;
plot(bins,lickhistoc_hit(14,:),'b');
plot(bins,lickhistr_hit(14,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,5);title('tone');hold on;
plot(bins,lickhistot_hit(10,:),'b');
plot(bins,lickhistr_hit(10,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,6);title('tone');hold on;
plot(bins,lickhistot_hit(15,:),'b');
plot(bins,lickhistr_hit(15,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,7);title('choice');hold on;
plot(bins,lickhistoc_hit(10,:),'b');
plot(bins,lickhistr_hit(10,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,8);title('choice');hold on;
plot(bins,lickhistoc_hit(15,:),'b');
plot(bins,lickhistr_hit(15,:),'k');
xline(2.8);xline(0.3);axis tight;
filetype='fig';
saveas(figAction,[subjlist{nbsubj} '_licks_after_light_Action' num2str(nbins) 'bins.' filetype]);


%% 29, 30, 32, 34 for sk204
nbsubj=7;
lickhistr_miss = lickhistrmiss{nbsubj};
lickhistr_cr = lickhistrcr{nbsubj};
lickhisto_miss= lickhistomiss{nbsubj};
lickhisto_cr = lickhistocr{nbsubj};
lickhistot_miss= lickhistotmiss{nbsubj};
lickhistot_cr=lickhistotcr{nbsubj};
lickhistoc_miss=lickhistocmiss{nbsubj};
lickhistoc_cr=lickhistoccr{nbsubj};


figlight=figure; title('sk204');
subplot(2,4,1);title('full trial');hold on;
plot(bins,lickhisto_miss(29,:),'b');
plot(bins,lickhistr_miss(29,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,2);title('full trial');hold on;
plot(bins,lickhisto_miss(32,:),'b');
plot(bins,lickhistr_miss(32,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,3);title('choice');hold on;
plot(bins,lickhistoc_miss(29,:),'b');
plot(bins,lickhistr_miss(29,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,4);title('choice');hold on;
plot(bins,lickhistoc_miss(32,:),'b');
plot(bins,lickhistr_miss(32,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,5);title('tone');hold on;
plot(bins,lickhistot_miss(30,:),'b');
plot(bins,lickhistr_miss(30,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,6);title('tone');hold on;
plot(bins,lickhistot_miss(34,:),'b');
plot(bins,lickhistr_miss(34,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,7);title('choice');hold on;
plot(bins,lickhistoc_miss(30,:),'b');
plot(bins,lickhistr_miss(30,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

subplot(2,4,8);title('choice');hold on;
plot(bins,lickhistoc_miss(34,:),'b');
plot(bins,lickhistr_miss(34,:),'k');
xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);
filetype='fig';
saveas(figlight,[subjlist{nbsubj} '_licks_after_light' num2str(nbins) 'bins.' filetype]);

%sk204
lickhistr_hit =lickhistrhit{nbsubj};
lickhistr_fa =lickhistrfa{nbsubj};
lickhisto_hit =lickhistohit{nbsubj};
lickhisto_fa= lickhistofa{nbsubj};
lickhistot_hit = lickhistothit{nbsubj};
lickhistot_fa = lickhistotfa{nbsubj};
lickhistoc_hit = lickhistochit{nbsubj};
lickhistoc_fa = lickhistocfa{nbsubj};
lickhistp_hit = lickhistphit{nbsubj};
lickhistp_fa = lickhistpfa{nbsubj};

figAction=figure; title('sk198');
subplot(2,4,1);title('full trial hit');hold on;
plot(bins,lickhisto_hit(29,:),'b');
plot(bins,lickhistr_hit(29,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,2);title('full trial hit');hold on;
plot(bins,lickhisto_hit(32,:),'b');
plot(bins,lickhistr_hit(32,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,3);title('choice');hold on;
plot(bins,lickhistoc_hit(29,:),'b');
plot(bins,lickhistr_hit(29,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,4);title('choice');hold on;
plot(bins,lickhistoc_hit(32,:),'b');
plot(bins,lickhistr_hit(32,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,5);title('tone');hold on;
plot(bins,lickhistot_hit(30,:),'b');
plot(bins,lickhistr_hit(30,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,6);title('tone');hold on;
plot(bins,lickhistot_hit(34,:),'b');
plot(bins,lickhistr_hit(34,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,7);title('choice');hold on;
plot(bins,lickhistoc_hit(30,:),'b');
plot(bins,lickhistr_hit(30,:),'k');
xline(2.8);xline(0.3);axis tight;

subplot(2,4,8);title('choice');hold on;
plot(bins,lickhistoc_hit(34,:),'b');
plot(bins,lickhistr_hit(34,:),'k');
xline(2.8);xline(0.3);axis tight;
filetype='fig';
saveas(figAction,[subjlist{nbsubj} '_licks_after_light_Action' num2str(nbins) 'bins.' filetype]);
