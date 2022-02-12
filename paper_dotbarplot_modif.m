function [out]=paper_dotbarplot_modif(x, xwt, colorm, sizedots, errorbarloc)


hold all
bar(x,nanmean(xwt), 0.4, 'Facecolor', 'none', 'Edgecolor', colorm, 'LineWidth',1.2)
plotSpread(xwt','xValues', x, 'distributionMarkers','o', 'distributionColors', colorm, 'spreadWidth',0.8,'sizedots', sizedots) %before spreadwithd was 1
set(gca,'fontsize',12)
if errorbarloc==1 %1 for SEM, 2 for STD and 0 for no error bar
    errorbar(x, nanmean(xwt), nanstd(xwt)/sqrt(length(xwt)), 'color', colorm, 'LineStyle','none','LineWidth',1.2);
elseif errorbarloc==2
    errorbar(x, nanmean(xwt), nanstd(xwt), 'color', colorm, 'LineStyle','none','LineWidth',1.2);
end
out=0;
end