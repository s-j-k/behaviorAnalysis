function [] = SingleAnimalPlots(datan, dataNew, animals)
%%% For individual animal
for i = 1:length(dataNew)
    clear xx_Hit yy_Hit xx_CR yy_CR xx_FA yy_FA LickLatencyFA LickLatencyHit
    figure(i)
    ColorMat={[96, 211, 148]./255; [242, 143, 59]./255};
    % Action rates, False alarm and Hit Rates for the Reinforced context
    subplot(3, 3, 1)
    plot(1:length(dataNew(i).reinfhit), (dataNew(i).reinfhit).*100,'Color',ColorMat{1}, 'LineWidth', 2), hold on;
    plot(1:length(dataNew(i).reinffa), (dataNew(i).reinffa).*100,'Color',ColorMat{2}, 'LineWidth', 2), hold on;
    title('Action Rates for Reinforced')
    %     legend
    ylim([0, 100])
    yticks([0 50 100]); yticklabels([0 50 100]);
    ylabel('Action Rate %')
    xlabel('Block of 20 trials')
    
    subplot(3,3, 2)
    ProbeHit = dataNew(i).probehit(~isnan(dataNew(i).probehit));
    ProbeFa = dataNew(i).probefa(~isnan(dataNew(i).probefa));
    plot(1:length(ProbeHit),ProbeHit.*100 ,'Color',ColorMat{1}, 'LineWidth', 2), hold on;
    plot(1:length(ProbeFa),ProbeFa.*100,'Color',ColorMat{2}, 'LineWidth', 2), hold on;
    title('Action Rates for Probe')
    ylim([0, 100])
    yticks([0 50 100]); yticklabels([0 50 100]);
    ylabel('Action Rate %')
    xlabel('Days')
    % %     legend
    
    subplot(3, 3, 3)
    dprimeProbe = dataNew(i).dprimeprobe(~isnan(dataNew(i).dprimeprobe));
    plot((1:length(dprimeProbe)), dprimeProbe, 'k','LineWidth', 2 ), hold on;
    plot((1:length(dataNew(i).dprimereinf)), dataNew(i).dprimereinf(:), 'b','LineWidth', 2), hold on;
    title('d prime')
    legend({'Probe', 'Reinforced'})
    ylabel('trials')
    xlabel('d prime')
    
    subplot(3, 3, 4)
    %     %      x = -2:0.5:2;
    %     w = 1;
    %     tr = 1;
    %     for m = 1:length(dataNew(i).HitLicks)
    %         for n = 1:length(dataNew(i).HitLicks{m})
    %             if ~isnan(dataNew(i).HitLicks{m}{n})
    %                 for p = 1:length(dataNew(i).HitLicks{m}{n})
    %                     xx_Hit(w)=  dataNew(i).HitLicks{m}{n}(p);
    %                     xx_Hit(w+1)= dataNew(i).HitLicks{m}{n}(p);
    %                     xx_Hit(w+2)= nan;
    %                     yy_Hit(w)= tr+0.5;
    %                     yy_Hit(w+1)= tr-0.5;
    %                     yy_Hit(w+2)= nan;
    %                     w = w+3;
    %                 end
    %                 tr = tr+1;
    %             else
    %             end
    %
    %         end
    %     end
    %     plot(xx_Hit, yy_Hit,'Color',ColorMat{1}, 'LineWidth', 1)
    %     title('Target Lick Rastor plot')
    %     xlabel('time(sec)')
    %     ylabel('#trials')
    edges = -2:0.05:2;
    %%% mean Lick latency trial wise
    Licks = [];
    LicksCount = [];
    for m = 1:length(dataNew(i).HitLicks)
        for  j= 1:length(dataNew(i).HitLicks{m})
            if ~isnan(dataNew(i).HitLicks{m}{j})
                Licks = (dataNew(i).HitLicks{m}{j});
                [xx yy] = histcounts(Licks, edges);
                LicksCount = [LicksCount; xx];
            else
            end
        end
    end
    meanLickRate = nanmean(LicksCount, 1);
    plot(edges(1:80), meanLickRate, 'b', 'LineWidth', 1.5), hold on;
    %     for k = 1:length(LicksCount)
    %     plot(edges(1:80), LicksCount(k,:)./(sum(LicksCount(k, :))), 'y'), hold on;
    %     end
    xlabel('time')
    ylabel('Mean Lick Counts')
    title('Hit licks in Reinforced')
    
    
    subplot(3, 3, 5)
%     q = 1;
%     tr = 1;
%     for m = 1:length(dataNew(i).FALicks)
%         for n = 1:length(dataNew(i).FALicks{m})
%             if ~isnan(dataNew(i).FALicks{m}{n})
%                 for p = 1:length(dataNew(i).FALicks{m}{n})
%                     xx_FA(q)=  dataNew(i).FALicks{m}{n}(p);
%                     xx_FA(q+1)= dataNew(i).FALicks{m}{n}(p);
%                     xx_FA(q+2)= nan;
%                     yy_FA(q)= tr+0.5;
%                     yy_FA(q+1)= tr-0.5;
%                     yy_FA(q+2)= nan;
%                     q = q+3;
%                 end
%                 tr = tr+1;
%             else
%             end
%         end
%     end
%     plot(xx_FA, yy_FA,'Color',ColorMat{2})
%     title('Foil Lick Rastor plot')
%     xlabel('time(sec)')
%     ylabel('#trials')
    edges = -2:0.05:2;
    %%% mean Lick latency trial wise
    Licks_FA = [];
    LicksCountFA = [];
    for m = 1:length(dataNew(i).FALicks)
        for  j= 1:length(dataNew(i).FALicks{m})
            if ~isnan(dataNew(i).FALicks{m}{j})
                Licks_FA = (dataNew(i).FALicks{m}{j});
                [xx yy] = histcounts(Licks_FA, edges);
                LicksCountFA = [LicksCountFA; xx];
            else
            end
        end
    end
    meanLickRateFA = nanmean(LicksCountFA, 1);
    plot(edges(1:80), meanLickRateFA, 'b', 'LineWidth', 1.5), hold on;
    %     for k = 1:length(LicksCount)
    %     plot(edges(1:80), LicksCount(k,:)./(sum(LicksCount(k, :))), 'y'), hold on;
    %     end
    xlabel('time')
    ylabel('Mean Lick Counts')
    title('FA licks in Reinforced')
    
    subplot(3, 3, 6)
    %     LickReinfCR = datan(i).lick(find(cell2mat(datan(i).outcome)==4));
%     w = 1;
%     tr = 1;
%     for m = 1:length(dataNew(i).CRLicks)
%         for n = 1:length(dataNew(i).CRLicks{m})
%             if ~isnan(dataNew(i).CRLicks{m}{n})
%                 for p = 1:length(dataNew(i).CRLicks{m}{n})
%                     xx_CR(w)= dataNew(i).CRLicks{m}{n}(p);
%                     xx_CR(w+1)= dataNew(i).CRLicks{m}{n}(p);
%                     xx_CR(w+2)= nan;
%                     yy_CR(w)= tr+0.5;
%                     yy_CR(w+1)= tr-0.5;
%                     yy_CR(w+2)= nan;
%                     w = w+3;
%                 end
%                 tr = tr+1;
%             else
%             end
%         end
%     end
%     %     figure(5)
%     plot(xx_CR, yy_CR,'k')
%     title('Correct reject Rastor plot')
%     ylabel('#trials/100')
%     xlabel('time(s)')
       edges = -2:0.05:2;
    %%% mean Lick latency trial wise
    Licks_CR = [];
    LicksCountCR = [];
    for m = 1:length(dataNew(i).CRLicks)
        for  j= 1:length(dataNew(i).CRLicks{m})
            if ~isnan(dataNew(i).CRLicks{m}{j})
                Licks_CR = (dataNew(i).CRLicks{m}{j});
                [xx yy] = histcounts(Licks_CR, edges);
                LicksCountCR = [LicksCountCR; xx];
            else
            end
        end
    end
    meanLickRateCR = nanmean(LicksCountCR, 1);
    plot(edges(1:80), meanLickRateCR, 'b', 'LineWidth', 1.5), hold on;
    %     for k = 1:length(LicksCount)
    %     plot(edges(1:80), LicksCount(k,:)./(sum(LicksCount(k, :))), 'y'), hold on;
    %     end
    xlabel('time')
    ylabel('Mean Lick Counts')
    title('Correct Rejects Licks')
      
    subplot(3, 3, 7)
    %%% finding the lick latency
    tr = 1;
    for k = 1:length(dataNew(i).HitLicks)
        for j = 1:length(dataNew(i).HitLicks{k})
            if ~isnan(dataNew(i).HitLicks{k}{j})
                idx = find(dataNew(i).HitLicks{k}{j}>0);
                LickLatencyHit(tr)= dataNew(i).HitLicks{k}{j}(idx(1));
                tr = tr+1;
            else
            end
        end
    end
    histogram(LickLatencyHit)
    title('Lick Latency for Hit trials')
    xlabel('Time(s)')
    ylabel('Frequency')
    
    
    subplot(3,3,8)
    tr = 1;
    for k = 1:length(dataNew(i).FALicks)
        for j = 1:length(dataNew(i).FALicks{k})
            if ~isnan(dataNew(i).FALicks{k}{j})
                idx = find(dataNew(i).FALicks{k}{j}>0);
                LickLatencyFA(tr)= dataNew(i).FALicks{k}{j}(idx(1));
                tr = tr+1;
            else
            end
        end
    end
    histogram(LickLatencyFA)
    title('Lick latency for false alarm trials')
    xlabel('Time(s)')
    ylabel('Frequency')
    

    
    sgtitle(animals(i))
end