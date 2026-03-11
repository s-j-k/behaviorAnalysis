function []=CatchData(rates,allDataTestsOnly,lickMatTestsOnly,reinfcolor,optocolor)

nbins = 50;
bins = linspace(-1,4,nbins); 

plotSingleAnimal=0;
if plotSingleAnimal==0
else
    clear uu;uu=2; %plot each animal
    for ww=1:length(allDataTestsOnly)
        animalRates=rates{ww,2}; % which animal are we plotting
        yyyFig=figure('Position', [100 100 800 200]);
        subplot(1,3,1);hold on;
        end1=size(animalRates,1);
        yyy=bar([animalRates(end1-uu,1) animalRates(end1-uu,5); ...
            animalRates(end1-uu,2) animalRates(end1-uu,6); 0 animalRates(end1-uu,10)]);
        % rates(i,:) = [rhit rfa phit pfa ohit ofa ...
        %             rfhit rffa cofffa conffa phit2 pfa2 ...
        %             othit otfa ochit ocfa];
        yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
        yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
        yyy(2).CData=[optocolor;optocolor;optocolor];box off; 
        xticklabels({'hit','','fa','','catch fa'}); legend('off','on');
        title([string(allDataTestsOnly{animalIdx(ww),1}) 'MGB Catch']);

        uu=1;subplot(1,3,2); hold on;
        end1=size(animalRates,1);
        yyy=bar([animalRates(end1-uu,1) animalRates(end1-uu,5); ...
            animalRates(end1-uu,2) animalRates(end1-uu,6); 0 animalRates(end1-uu,10)]);
        yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
        yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
        yyy(2).CData=[optocolor;optocolor;optocolor]; box off; 
        xticklabels({'hit','fa','catch fa'}); legend('off','on');
        title('IC Catch');

        subplot(1,3,3);
        end1=size(animalRates,1);
        yyy=bar([animalRates(end1,1) animalRates(end1,5); ...
            animalRates(end1,2) animalRates(end1,6); 0 animalRates(end1,10)]);
        yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
        yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
        yyy(2).CData=[optocolor;optocolor;optocolor]; box off; 
        xticklabels({'hit','fa','catch fa'}); legend('off','on');
        title('LOB Catch');       
        animalName=allDataTestsOnly{animalIdx(ww),1};
            saveas(yyyFig,[char(allDataTestsOnly{animalIdx(ww),1}) '-CatchBar.fig']);
            saveas(yyyFig,[char(allDataTestsOnly{animalIdx(ww),1}) '-CatchBar.png']);
            close(yyyFig);

    end
end

% average of all animals
uu=2;
for ww=2:length(rates)
    animalRates=rates{ww,2}; 
    end1=size(animalRates,1);
    rhit(ww)=animalRates(end1-uu,1);
    rfa(ww)=animalRates(end1-uu,2);
    coffa(ww)=animalRates(end1-uu,10);
    % rates(i,:) = [rhit rfa phit pfa ohit ofa ...
    %             rfhit rffa cofffa conffa phit2 pfa2 ...
    %             othit otfa ochit ocfa];
end
yyyFig=figure('Position', [100 100 800 200]);
hold on;
subplot(1,2,1);
rfa=rfa(2:7);coffa=coffa(2:7);
rfa=rfa(~isnan(coffa));
coffa=coffa(~isnan(coffa));
yy=bar([nanmean(rfa) nanmean(coffa)]); hold on;
scatter(repmat(yy(1).XEndPoints(1),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(yy(1).XEndPoints(2),size(coffa,2),1), ...
    coffa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(mean(rfa),mean(coffa));
sigstar({[1,2]}, p)
yy(1).FaceColor='flat';
yy(1).CData=[reinfcolor;optocolor];ylim([0 1]);
xticklabels({'fa','catch fa'}); legend('off','on');
title(['Averaged MGB Catch']);
ylabel('Action rate');


catFa=[];catCr=[];rMiss=[];rFa=[];rCr=[];
for uu=2:size(lickMatTestsOnly,1)
    catch_fa=lickMatTestsOnly{uu,3};
    lickhistr_miss = lickMatTestsOnly{uu,14};
    lickhistr_fa =lickMatTestsOnly{uu,5};
    lickhistr_cr = lickMatTestsOnly{uu,15};
    catFa=vertcat(catFa,catch_fa);
    rMiss=vertcat(rMiss,lickhistr_miss);
    rFa=vertcat(rFa,lickhistr_fa);
    rCr=vertcat(rCr,lickhistr_cr);
end

% lick histograms
subplot(1,2,2);
shadedErrorBar(bins,nanmean(rFa,1),std(rFa,1,"omitnan"),'k');hold on;
shadedErrorBar(bins,nanmean(catFa,1),std(catFa,1,"omitnan"),'b');
    [h,pHit,ci,stats] = ttest2(nanmean(catFa,1),nanmean(rFa,1));
    text(min(xlim), min(ylim),num2str(pHit));
    xline(2.8);xline(0.3);axis tight;ylim([0 1.5]);xlim([0 4]);
    box off; title('Catch FA');
ylim([0 1]);ylabel('Distribution of licks');xlabel('Time in trial (s)');
yyyFig.Position(3:4)=[525 200];
saveas(gcf,['By Animal_T_MGB_catchHist_Opto']);
saveas(gcf,['By Animal_T_MGB_catchHist_Opto.png']);
saveas(gcf,['By Animal_T_MGB_catchHist_Opto.pdf']);

    
end