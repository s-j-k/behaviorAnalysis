t = 1;
ColorMat={[96, 211, 148]./255; [242, 143, 59]./255};
 nrow = length(datan);
for k = 1:length(datan) %% animal number
%for k = 1 %
    %%
    clear idxhitIndDay idxFAIndDay idxFAIndDayDiff FAIndDayIdx IndDayFA idxhitIndTrial idxFAIndTrial FAIndTrialIdx
    subplot(nrow, 3, t)
    plot(1:length(DataDay(k).reinfhit),(DataDay(k).reinfhit).*100 , 'Color', ColorMat{1}, 'linewidth', 2), hold on
    plot(1:length(DataDay(k).reinffa),(DataDay(k).reinffa).*100 , 'Color', ColorMat{2}, 'linewidth', 2),
    ylim([0, 100])
    
    idxStartDay = find((DataDay(k).reinfhit).*100>=80 & ((DataDay(k).reinffa).*100>=70));
    plot(idxStartDay(1)*ones(1, 100), 1:100, '--k', 'linewidth', 2), hold on
    DayStart = idxStartDay(1);
    
    
    idxhitIndDay = find((DataDay(k).reinfhit).*100>80);
    idxFAIndDay = find(((DataDay(k).reinffa).*100)<20);
    idxFAIndDayDiff = diff(idxFAIndDay);
    FAIndDayIdx = [];
    for i = 1:length(idxFAIndDayDiff)
        if idxFAIndDayDiff(i)==1
            if length(idxFAIndDayDiff)==1
                FAIndDayIdx = idxFAIndDay(1);
            elseif idxFAIndDayDiff(i+1)==idxFAIndDayDiff(i)==1
                FAIndDayIdx = idxFAIndDay(i);
                break
            else
            end
        else
            continue
        end
    end
    if ~isnan(FAIndDayIdx)
    IndDayFA =(FAIndDayIdx);
    plot(IndDayFA*ones(1, 100), 1:100, '--k', 'linewidth', 2), hold on
    %text(IndDayFA+1, 50, ['Day:' num2str(IndDayFA)], 'FontSize', 14, 'FontWeight', 'Bold');
    text(IndDayFA+1, 50, ['Day:' num2str(IndDayFA-DayStart)], 'FontSize', 14, 'FontWeight', 'Bold');
    else
        IndDayFA = NaN;
    end
    
    ylim([0, 100])
    % ErrorBand(meanReinfHit, sErrorMeanHit, x, ColorMat{1});
    % ErrorBand(meanReinfFA, sErrorMeanFA,x,  ColorMat{2});
    % title('Action Rates for Reinforced')
    bx = gcf;
    bx.Position(3:4) = [590 270];
    % set(gcf, 'position', [0.13, 0.11, 0.3, 0.6])
    ylim([0, 100])
    % xlim([0, (length(dataNew(1).reinfhit)*70)])
    yticks([0 50 100]); yticklabels([0 50 100]);
    % ylabel('Action Rate')
    % xlabel('Trials')
    ax = gca;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    ax.XLabel.FontSize = 14;
    ax.YLabel.FontSize = 14;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    
    %%
    subplot(nrow, 3, t+1)
    x = ([1 ((1:(length(dataNew(k).reinfhit)-1)).*70+1)]);
    % plot(1:length(meanReinfHit), meanReinfHit(:), 'Color', ColorMat{1}, 'linewidth', 2), hold on;
    % plot(1:length(meanReinfFA), meanReinfFA(:), 'Color', ColorMat{2}, 'linewidth', 2), hold on;
    plot(x, dataNew(k).reinfhit.*100, 'Color', ColorMat{1}, 'linewidth', 2), hold on;
    plot(x, dataNew(k).reinffa.*100, 'Color', ColorMat{2}, 'linewidth', 2), hold on;
    
    idxStartTrial = find((dataNew(k).reinfhit).*100>=80 & (((dataNew(k).reinffa).*100>=70)));
    
%     if idxStartTrial(1)==1
%     TrialStart = (idxStartTrial(1));
%     else
%         TrialStart = idxStartTrial(1)*70+1;
%     end
%     
%     plot(TrialStart*ones(1, 100), 1:100, '--k', 'linewidth', 2), hold on
%     
%     idxhitIndTrial = find((dataNew(k).reinfhit).*100>80);
%     idxFAIndTrial = find(((dataNew(k).reinffa).*100)<20);
%     idxFAIndTrial = idxFAIndTrial(idxFAIndTrial>10);
%     idxFAIndTrialDiff = diff(idxFAIndTrial);
    %     idxFAIndTrialDiff = idxFAIndTrialDiff(5:end);
%     FAIndTrialIdx = [];
%     for i = 1:length(idxFAIndTrialDiff)
%         if idxFAIndTrialDiff(i)==1
%             if idxFAIndTrialDiff(i+1)==idxFAIndTrialDiff(i)==1
%                 FAIndTrialIdx = idxFAIndTrial(i);
%                 break
%             else
%             end
%         else
%             continue
%         end
%     end
    if ~isnan(FAIndTrialIdx)
    IndTrialFA =(FAIndTrialIdx)*70+1;
    plot(IndTrialFA*ones(1, 100), 1:100, '--k', 'linewidth', 2), hold on
    %text(IndTrialFA+1, 50, ['Trial:' num2str(IndTrialFA)], 'FontSize', 14, 'FontWeight', 'Bold');
    text(IndTrialFA+1, 50, ['Trials:' num2str(IndTrialFA- TrialStart)], 'FontSize', 14, 'FontWeight', 'Bold');
    else
        IndTrialFA = NaN;
    end
    
    ylim([0, 100])
    % ErrorBand(meanReinfHit, sErrorMeanHit, x, ColorMat{1});
    % ErrorBand(meanReinfFA, sErrorMeanFA,x,  ColorMat{2});
    % title('Action Rates for Reinforced')
    bx = gcf;
    bx.Position(3:4) = [590 270];
    % set(gcf, 'position', [0.13, 0.11, 0.3, 0.6])
    ylim([0, 100])
    xlim([0, (length(dataNew(1).reinfhit)*70)])
    yticks([0 50 100]); yticklabels([0 50 100]);
    % ylabel('Action Rate')
    % xlabel('Trials')
    ax = gca;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    ax.XLabel.FontSize = 14;
    ax.YLabel.FontSize = 14;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    
    
    %%
    subplot(nrow, 3, t+2)
    plot(1:length(DataDay(k).probehit), DataDay(k).probehit.*100,  'Color',ColorMat{1}, 'linewidth', 2), hold on
    plot(1:length((DataDay(k).probefa)), DataDay(k).probefa.*100, 'Color',ColorMat{2}, 'linewidth', 2), hold on
    bx = gcf;
    bx.Position(3:4) = [450 200];
    ylim([0, 100])
    %xlim([0, (length(dataNew(1).reinfhit)*70)])
    yticks([0 50 100]); yticklabels([0 50 100]);
    % ylabel('Action Rate')
    % xlabel('Trials')
    ax = gca;
    ax.XLabel.FontWeight = 'bold';
    ax.YLabel.FontWeight = 'bold';
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontWeight = 'bold';
    ax.XLabel.FontSize = 14;
    ax.YLabel.FontSize = 14;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
    
    
    t = t+3;
    %% Extracting data to an excel file 
    %fileId1 = fopen('Animals_HitFAData.xls', 'a');
    cd 'K:/Behavior/'
    fileId2 = fopen('Animals_HitFAData.txt', 'a+');
    fprintf(fileId2, '%s, %i, %i, %i, %i, %i, %i\n', animals{k}, DayStart, IndDayFA, IndDayFA-DayStart, TrialStart, IndTrialFA, IndTrialFA-TrialStart)
    fclose(fileId2)
end