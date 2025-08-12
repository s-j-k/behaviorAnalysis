function plotlicks_after_light(lickMatTestsOnly)

% for sk198 39+42 - full and choice, 40+44 - tone and choice
nbins = 50;
bins = linspace(-1,4,nbins); 
rfa=[];ofa=[];cfa=[];tfa=[];
rhit=[];ohit=[];chit=[];thit=[];
rmiss=[];omiss=[];cmiss=[];tmiss=[];
rcr=[];ocr=[];ccr=[];tcr=[];
for uu=2:size(lickMatTestsOnly,1)
    lickhistr_fa =lickMatTestsOnly{uu,5};
    lickhisto_hit =lickMatTestsOnly{uu,6};
    lickhisto_fa= lickMatTestsOnly{uu,7};
    lickhistot_hit = lickMatTestsOnly{uu,8};
    lickhistot_fa = lickMatTestsOnly{uu,9};
    lickhistoc_hit = lickMatTestsOnly{uu,10};
    lickhistoc_fa = lickMatTestsOnly{uu,11};    
    lickhistr_miss = lickMatTestsOnly{uu,14};
    lickhistr_cr = lickMatTestsOnly{uu,15};
    lickhisto_miss= lickMatTestsOnly{uu,16};
    lickhisto_cr = lickMatTestsOnly{uu,17};
    lickhistot_miss= lickMatTestsOnly{uu,18};
    lickhistot_cr=lickMatTestsOnly{uu,19};
    lickhistoc_miss=lickMatTestsOnly{uu,20};
    lickhistoc_cr=lickMatTestsOnly{uu,21};
    lickhistr_hit =lickMatTestsOnly{uu,4};

    mousenum=lickMatTestsOnly{uu,1};
    figlight=figure; title(mousenum);
    
    
    subplot(4,3,1);title('full trial hit');hold on;
    plot(bins,nanmean(lickhisto_hit,1),'b');
    plot(bins,nanmean(lickhistr_hit,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,2);title('choice hit');hold on;
    plot(bins,nanmean(lickhistoc_hit,1),'b');
    plot(bins,nanmean(lickhistr_hit,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,3);title('tone hit');hold on;
    plot(bins,nanmean(lickhistot_hit,1),'b');
    plot(bins,nanmean(lickhistr_hit,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,4);title('full trial, miss,');hold on;
    plot(bins,nanmean(lickhisto_miss,1),'b');
    plot(bins,nanmean(lickhistr_miss,1),'k');
    xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

    subplot(4,3,5);title('choice, miss');hold on;
    plot(bins,nanmean(lickhistoc_miss,1),'b');
    plot(bins,nanmean(lickhistr_miss,1),'k');
    xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

    subplot(4,3,6);title('tone, miss');hold on;
    plot(bins,nanmean(lickhistot_miss,1),'b');
    plot(bins,nanmean(lickhistr_miss,1),'k');
    xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

    subplot(4,3,7);title('full trial cr');hold on;
    plot(bins,nanmean(lickhisto_cr,1),'b');
    plot(bins,nanmean(lickhistr_cr,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,8);title('choice cr');hold on;
    plot(bins,nanmean(lickhistoc_cr,1),'b');
    plot(bins,nanmean(lickhistr_cr,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,9);title('tone cr');hold on;
    plot(bins,nanmean(lickhistot_cr,1),'b');
    plot(bins,nanmean(lickhistr_cr,1),'k');
    xline(2.8);xline(0.3);axis tight;
    
    subplot(4,3,10);title('full trial fa');hold on;
    plot(bins,nanmean(lickhisto_fa,1),'b');
    plot(bins,nanmean(lickhistr_fa,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,11);title('choice fa');hold on;
    plot(bins,nanmean(lickhistoc_fa,1),'b');
    plot(bins,nanmean(lickhistr_fa,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,12);title('tone fa');hold on;
    plot(bins,nanmean(lickhistot_fa,1),'b');
    plot(bins,nanmean(lickhistr_fa,1),'k');
    xline(2.8);xline(0.3);axis tight;

    filetype='fig';
    saveas(figlight,[mousenum '_licks_after_light' num2str(nbins) 'bins.' filetype]);
    filetype='pdf';
    saveas(figlight,[mousenum '_licks_after_light' num2str(nbins) 'bins.' filetype]);
    
    rfa=vertcat(rfa,lickhistr_fa);
    ohit= vertcat(ohit,lickhisto_hit);
    ofa=vertcat(ofa,lickhisto_fa);
    cfa=vertcat(cfa,lickhistoc_fa);
    tfa=vertcat(tfa,lickhistot_fa);
    rhit=vertcat(rhit,lickhistr_hit);
    chit=vertcat(chit,lickhistoc_hit);
    thit=vertcat(thit,lickhistot_hit);
    rmiss=vertcat(rmiss,lickhistr_miss);
    omiss=vertcat(omiss,lickhisto_miss);
    cmiss=vertcat(cmiss,lickhistoc_miss);tmiss=vertcat(tmiss,lickhistot_miss);
    rcr=vertcat(rcr,lickhistr_cr);ocr=vertcat(ocr,lickhisto_cr);
    ccr=vertcat(ccr,lickhistoc_cr);tcr=vertcat(tcr,lickhistot_cr);
    
end

    figlightAvg=figure;
    
    subplot(4,3,1);title('full trial hit');hold on;
    plot(bins,nanmean(ohit,1),'b');
    plot(bins,nanmean(rhit,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,2);title('choice hit');hold on;
    plot(bins,nanmean(chit,1),'b');
    plot(bins,nanmean(rhit,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,3);title('tone hit');hold on;
    plot(bins,nanmean(thit,1),'b');
    plot(bins,nanmean(rhit,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,4);title('full trial, miss,');hold on;
    plot(bins,nanmean(omiss,1),'b');
    plot(bins,nanmean(rmiss,1),'k');
    xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

    subplot(4,3,5);title('choice, miss');hold on;
    plot(bins,nanmean(cmiss,1),'b');
    plot(bins,nanmean(rmiss,1),'k');
    xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

    subplot(4,3,6);title('tone, miss');hold on;
    plot(bins,nanmean(tmiss,1),'b');
    plot(bins,nanmean(rmiss,1),'k');
    xline(2.8);xline(0.3);axis tight;ylim([0 0.3]);

    subplot(4,3,7);title('full trial cr');hold on;
    plot(bins,nanmean(ocr,1),'b');
    plot(bins,nanmean(rcr,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,8);title('choice cr');hold on;
    plot(bins,nanmean(ccr,1),'b');
    plot(bins,nanmean(rcr,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,9);title('tone cr');hold on;
    plot(bins,nanmean(tcr,1),'b');
    plot(bins,nanmean(rcr,1),'k');
    xline(2.8);xline(0.3);axis tight;
    
    subplot(4,3,10);title('full trial fa');hold on;
    plot(bins,nanmean(ofa,1),'b');
    plot(bins,nanmean(rfa,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,11);title('choice fa');hold on;
    plot(bins,nanmean(cfa,1),'b');
    plot(bins,nanmean(rfa,1),'k');
    xline(2.8);xline(0.3);axis tight;

    subplot(4,3,12);title('tone fa');hold on;
    plot(bins,nanmean(tfa,1),'b');
    plot(bins,nanmean(rfa,1),'k');
    xline(2.8);xline(0.3);axis tight;

    filetype='fig';
    saveas(figlightAvg,['AVG_allmice_licks_after_light' num2str(nbins) 'bins.' filetype]);
    filetype='pdf';
    saveas(figlightAvg,['AVG_allmice_licks_after_light' num2str(nbins) 'bins.' filetype]);
    
    

end
