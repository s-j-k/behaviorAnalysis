function [out]=paper_dotbarplot_new(xwt, colorm, xko, colorml, titlestr, want_outlier)

if want_outlier==1
    
    normthrwt=kstest(xwt(~Isoutlier(xwt)));
    normthrko=kstest(xko(~Isoutlier(xko)));
    if normthrwt==0 && normthrko==0
        [h pval]=ttest2(xwt(~Isoutlier(xwt)),xko(~Isoutlier(xko))); %if data are normally distributed (0)
    else
        [pval h]=ranksum(xwt(~Isoutlier(xwt)),xko(~Isoutlier(xko))); %if data are not normally distributed (1)
    end
    
    hold all
    bar(1,nanmean(xwt(~Isoutlier(xwt))), 0.4, 'Facecolor', 'none', 'Edgecolor', colorm)
    bar(1.5,nanmean(xko(~Isoutlier(xko))), 0.4, 'Facecolor', 'none', 'Edgecolor', colorml)
    if ~isempty(xwt(Isoutlier(xwt)))
        plot(1,xwt(Isoutlier(xwt)), 'xk', 'color', colorm)
    end
    if ~isempty(xko(Isoutlier(xko)))
        plot(1.5,xko(Isoutlier(xko)), 'x', 'color', colorml)
    end
    plotSpread(xwt(~Isoutlier(xwt))','xValues', 1, 'distributionMarkers','.', 'distributionColors', colorm, 'spreadWidth',1, 'showMM', 4)
    plotSpread(xko(~Isoutlier(xko))','xValues', 1.5, 'distributionMarkers','.', 'distributionColors',colorml, 'spreadWidth',1, 'showMM', 4)
    set(gca,'fontsize',12)
    xlim([0.6 1.9])
    title([titlestr ' p=' num2str(pval)])
    set(gca, 'xtick', [1 1.5], 'xticklabel',{[ ]})
    out=0;
    
elseif want_outlier==0
    
    normthrwt=kstest(xwt);
    normthrko=kstest(xko);
    if normthrwt==0 && normthrko==0
        [h pval]=ttest2(xwt,xko); %if data are normally distributed (0)
    else
        [pval h]=ranksum(xwt,xko); %if data are not normally distributed (1)
    end
    
    hold all
    bar(1,nanmean(xwt), 0.4, 'Facecolor', 'none', 'Edgecolor', colorm)
    bar(1.5,nanmean(xko), 0.4, 'Facecolor', 'none', 'Edgecolor', colorml)
%     if ~isempty(xwt)
%         plot(1,xwt, 'xk')
%     end
%     if ~isempty(xko)
%         plot(1.5,xko, 'x', 'color', colorml)
%     end
    plotSpread(xwt','xValues', 1, 'distributionMarkers','.', 'distributionColors', colorm, 'spreadWidth',0.3, 'showMM', 4) %before spreadwithd was 1
    plotSpread(xko','xValues', 1.5, 'distributionMarkers','.', 'distributionColors',colorml, 'spreadWidth',0.3, 'showMM', 4)
    set(gca,'fontsize',12)
    xlim([0.6 1.9])
    title([titlestr ' p=' num2str(pval)])
    set(gca, 'xtick', [1 1.5], 'xticklabel',{[ ]})
    out=0;
end
end