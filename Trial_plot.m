% individual mouse summary plot

SingleAnimalPlots3(datan, dataNew, animals, DataDay);




%% Big summary plot

ColorMat={[96, 211, 148]./255; [242, 143, 59]./255};
j = 1;
kk=length(animals);
ff=4;

for i = 1:length(DataDay)
    
    if strcmp(datan(i).exp, 'Control')==1
        subplot(kk,ff, j)
        plot(1:length(DataDay(i).reinfhit), DataDay(i).reinfhit.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
        plot(1:length((DataDay(i).reinffa)), DataDay(i).reinffa.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        xlabel('Days')
        ylabel('Action Rate')
        ylim([50 100])
        ax = gca;
        
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';
        
        subplot(kk,ff,j+1)
        plot(1:length(dataNew(i).reinfhit), dataNew(i).reinfhit.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
        plot(1:length((dataNew(i).reinffa)), dataNew(i).reinffa.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        xlabel('Blocks of 20 trials')
        ylim([50 100])
        ax = gca;
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';
        
        subplot(kk,ff,j+2)
        m = round( length(dataNew(i).reinfhit).*(1/2));
        plot(1:m, dataNew(i).reinfhit(1:m).*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
        plot(1:m, dataNew(i).reinffa(1:m).*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        xlabel('Blocks of 20 trials')
        ylim([50 100])
        ax = gca;
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';
        
        subplot(kk,ff,j+3)
        plot(1:length(DataDay(i).probehit), DataDay(i).probehit.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
        plot(1:length((DataDay(i).probefa)), DataDay(i).probefa.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        xlabel('Days')
        ylim([50 100])      
        ax = gca;
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';
    else
        subplot(kk,ff, j)
        plot(1:length(DataDay(i).reinfhit), DataDay(i).reinfhit.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
        plot(1:length((DataDay(i).reinffa)), DataDay(i).reinffa.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        xlabel('Days')
        ylabel('Action Rate Reinf.')
        ylim([50 100])       
        ax = gca;
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';
        
        
        subplot(kk,ff,j+1)
%         plot(1:length(dataNew(i).reinfhit), dataNew(i).reinfhit.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
%         plot(1:length((dataNew(i).reinffa)), dataNew(i).reinffa.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
%         xlabel('Blocks of 20 trials')
%         ylim([50 100])
%         ax = gca;
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';
        
        
%         subplot(kk,ff,j+2)
                m = round( length(dataNew(i).reinfhit).*(1/2));
                plot(1:m, dataNew(i).reinfhit(1:m),  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
                plot(1:m, dataNew(i).reinffa(1:m), 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        plot(1:length(dataNew(i).reinfhit), dataNew(i).reinfhit.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
        plot(1:length((dataNew(i).reinffa)), dataNew(i).reinffa.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        xlabel('Blocks of 20 trials')
        legend('R-Hit','R-FA')
        ylabel('Action Rate Reinf.')
%         ylim([50 100])       
%         ax = gca;
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';
        

        subplot(kk,ff,j+2)
%         subplot(kk,ff,j+3)
        plot(1:length(DataDay(i).probehit), DataDay(i).probehit.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
        plot(1:length((DataDay(i).probefa)), DataDay(i).probefa.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        xlabel('Days')
        legend('Probe Hit','Probe FA')
        ylabel('Action Rate Probe')
%         ylim([50 100])        
%         ax = gca;
%         ax.XLabel.FontWeight = 'bold';
%         ax.YLabel.FontWeight = 'bold';
%         ax.XAxis.FontWeight = 'bold';
%         ax.YAxis.FontWeight = 'bold';



        subplot(kk,ff,j+3)
        plot(1:length((dataNew(i).dprimereinf)), dataNew(i).dprimereinf.*100, 'Color',ColorMat{2}, 'linewidth', 1.5), hold on
        plot(1:length(DataDay(i).dprimeprobe), DataDay(i).dprimeprobe.*100,  'Color',ColorMat{1}, 'linewidth', 1.5), hold on
%         plot(1:length(dprime1),dprime1(1,:));
%         plot(1:length(dprime1),dprime1(2,:))
        legend('Reinforced','Probe')
        xlabel('days of training');ylabel('d`');
%         title(data(1,ii).mouse)


    end
    j = j+4;
end