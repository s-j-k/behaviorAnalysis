%For bpod data lick raster plot plotting for all training days of each
%mouse. Data needs to be preprocessed. 
%S.M. 12/03/2019
% clear all
% addpath('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Behavior\CitricAcid\')
% addpath('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Matlab_functions\Outside');

% load('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Behavior\CitricAcid\Analysis\data_GNGNormal_20200925_160827.mat')
% load('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Behavior\CitricAcid\Analysis\data_GNGNormal_20201115_merged.mat')

% load('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Behavior\CitricAcid\Analysis\data_GNGNormal_20210628_3cohorts.mat')

%To plot rasterplots of licks per animal per day
clear Licks LickMat Probetrials Reinforcedtrials 
% load('/Users/shar/Documents/APP_project/Analysis/APPdata_all_FINAL_corr2.mat')
wt=[];
mut=[];
ext=[];
for ii=1:length(data)
    Animals(ii,:)=data(ii).mouse;
    
    if data(ii).genotype==1
        wt=[wt ii];
%     elseif data(ii).genotype==2
%         mut=[mut ii];
%     elseif data(ii).genotype==3
%         ext=[ext ii];
    end 
end


for ii=1:length(data)
    mouse{ii}=data(ii).mouse;
    xmt=data(ii).lickstime; 
    Licks=xmt;
    for di=1:length(data(ii).expday)
%     for di=1
        ProbeTrials=find(data(ii).trialhist{1,di}>=3);
        ProbeTTargets=find(data(ii).trialhist{1,di}==3);
        ProbeTFoils=find(data(ii).trialhist{1,di}==4);
        ReinforcedTrials=find(data(ii).trialhist{1,di}<=2);
        Licksday=Licks{1,di};
        Targets=find(data(ii).trialhist{1,di}==1);
        Foils=find(data(ii).trialhist{1,di}==2);
                
        %for PSTH
        clear LickMat
        Tx=ReinforcedTrials;
        binsize=0.05; %in seconds. 0.05 is the one where most peaks are detected for reliable licking frequency
        bins=[-1:binsize:5+binsize];
        timevect=[];
        LickMat=NaN(length(Tx),length(bins)-2);
        for ti=1:length(Tx)
            for bi=1:length(bins)-2
                LickMat(Tx(ti),bi)=length(find(Licksday{1,Tx(ti)}>=bins(bi) & Licksday{1,Tx(ti)}<bins(bi+1)));
            end
        end
        
        LickDayMat_Targets{ii}(di,:)=sum(LickMat(Targets,:));
        if size(LickMat(Foils,:),1)>1
            LickDayMat_Foils{ii}(di,:)=sum(LickMat(Foils,:));
        elseif size(LickMat(Foils,:),1)==1
            LickDayMat_Foils{ii}(di,:)=LickMat(Foils,:);
        end
        HitsR{ii}(di,:)=data(ii).reinforced{1,di}(1)/length(Targets);
        FalseAlarmsR{ii}(di, :)=data(ii).reinforced{1,di}(4)/length(Foils);
        HitsP{ii}(di,:)=data(ii).probe{1,di}(1)/length(ProbeTTargets);
        FalseAlarmsP{ii}(di, :)=data(ii).probe{1,di}(4)/length(ProbeTFoils);
        
        dprimeR{ii}(di,:)=data(ii).dprime{1,di}(1);
        dprimeP{ii}(di,:)=data(ii).dprime{1,di}(2);
        
%         %to use only higher dprimes during probe 
%         if isnan(dprimeP{ii}(di,:))
%             dprimePtouse{ii}(di,:)=NaN;
%         elseif ~isnan(dprimeP{ii}(di,:))
%             if sum(isnan(dprimePtouse{ii}))==length(dprimePtouse{ii}) %if this is the first numeric d prime to be added
%                 dprimePtouse{ii}(di,:)=dprimeP{ii}(di,:);
%             elseif dprimeP{ii}(di,:)>=max(dprimePtouse{ii}) %if other max dprimes exist
%                 dprimePtouse{ii}(di,:)=dprimeP{ii}(di,:);
%             else
%                 dprimePtouse{ii}(di,:)=max(dprimePtouse{ii});
%             end
%         end
        
    end
end
%% To plot foil and target licks per day per animal
colorsmat2=[[155, 204, 85]./255; [228, 87, 46]./255; [179, 201, 174]./255; [239, 149, 137]./255];
% colorsmat2=[[179, 201, 174]./255; [239, 149, 137]./255];

for ii=1:length(data)
    xmt=data(ii).lickstime; %data(ii).lickstime; lickstimefirst
    Licks=xmt;
    for di=2 %:length(data(ii).expday)
        ProbeTrials=find(data(ii).trialhist{1,di}>=3);
        ReinforcedTrials=find(data(ii).trialhist{1,di}<=2);
        Licksday=Licks{1,di};
        Targets=find(data(ii).trialhist{1,di}==1);
        Foils=find(data(ii).trialhist{1,di}==2);
        
        %for PSTH
        Tx=ReinforcedTrials;
        binsize=0.05; %in seconds. 0.05 is the one where most peaks are detected for reliable licking frequency
        bins=[0:binsize:5+binsize];
        LickMat=NaN(length(Tx),length(bins)-2);
        for ti=1:length(Tx)
            for bi=1:length(bins)-2
                LickMat(Tx(ti),bi)=length(find(Licksday{1,Tx(ti)}>=bins(bi) & Licksday{1,Tx(ti)}<bins(bi+1)));
            end
        end
        
        figure
        subplot(2,1,1)
        hold all
        for ti=ReinforcedTrials
            if ismember(ti,Targets)
                colort=colorsmat2(1,:);
            elseif ismember(ti,Foils)
                colort=colorsmat2(2,:);
            end
            plot(Licks{1,di}{1,ti}, ti*ones(1,length(Licks{1,di}{1,ti})), '.', 'color', [colort])
        end
        %         if ~isempty(ProbeTrials)
        %             for ti=ProbeTrials
        %                 plot(Licks{1,di}{1,ti}, ti*ones(1,length(Licks{1,di}{1,ti})), '.', 'color', [98, 121, 184]]./255)
        %             end
        %         end
        line([0 0], [0 350], 'color', 'r');
        ylabel('Trial #')
        ylabel('Time')
        xlim([-0.5 4]); ylim([0 350])
        if ismember(ii,wt)
%             title(['\color{Blue}' char(Animals(ii,:)) ' / ' num2str(di)])
            title(['\color{Blue}' char(Animals(1,ii)) ' / ' num2str(di)])
%         elseif ismember(ii,ext)
%             title(['\color{Black}' char(Animals(ii,:)) ' / ' num2str(di)])
%         elseif ismember(ii,mut)
%             title(['\color{Yellow}' char(Animals(ii,:)) ' / ' num2str(di)])
        end
        
        subplot(2,1,2)
        hold all
        plot(0:length(nansum(LickMat(Targets,:)))-1,nansum(LickMat(Targets,:)), 'color', colorsmat2(1,:), 'LineWidth', 2)
        plot(0:length(nansum(LickMat(Foils,:)))-1,nansum(LickMat(Foils,:)), 'color', colorsmat2(2,:), 'LineWidth', 2)
        line([0 0], [0 max(nansum(LickMat(Targets,:)))+10], 'color', 'r');
        ylim([0 max(nansum(LickMat(Targets,:)))+10])
        xlim([-10 80])
        ylabel('Lick sum (50 ms bins)')
        xlabel('Bins')
        dpr=num2str(data(ii).dprime{1,di}(1));
        dpp=num2str(data(ii).dprime{1,di}(2));
        text(40, max(nansum(LickMat(Targets,:)))+5, ['d-reinf= ' dpr]);
        text(60, max(nansum(LickMat(Targets,:)))+5, ['d-prob= ' dpp]);
    end
end
%% Heatmaps of all training days and plots of performance
for ii=1:length(data)
    h=figure(ii);
    h.Position=[50 380 1300 300];
    h1=subplot(1,5,1);
    h1.Position=[0.05 0.11 0.24 0.81];
    hold on
    imagesc(LickDayMat_Targets{ii})
    caxis([0 max(max(LickDayMat_Targets{ii}))])
    ch1=colorbar; ch1.Position=[0.295 0.11 0.008 0.81];
    colormap(hot);
    ylim([0.5 size(LickDayMat_Targets{ii}, 1)+0.5])
    if ismember(ii,wt)
        ylabel(['\color{Blue}' char(mouse{ii}) ' / Day #' ])
    elseif ismember(ii,mut)
        ylabel([char(mouse{ii}) ' / Day #' ])
    end
    xlabel(['Time (s)'])
    title('Target trials')
    %line([find(bins==0)-0.5 find(bins==0)-0.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([20.5 20.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([30.5 30.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([40.5 40.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([50.5 50.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    xlimtarg=[0.5:4:120.5];
    xticks(xlimtarg)
    xticklabels(bins(1:4:end))
    xlim([16.5 60.5])
    hold all
    
    h2=subplot(1,5,2);
    h2.Position=[0.34 0.11 0.24 0.81];
    hold all
    imagesc(LickDayMat_Foils{ii})
    caxis([0 max(max(LickDayMat_Targets{ii}))])
    ch2=colorbar; ch2.Position=[0.585 0.11 0.008 0.81];
    colormap(hot);
    ylim([0.5 size(LickDayMat_Foils{ii}, 1)+0.5])
    %line([find(bins==0)-0.5 find(bins==0)-0.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([20.5 20.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([30.5 30.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([40.5 40.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    line([50.5 50.5], [0 size(LickDayMat_Targets{ii}, 1)+0.5], 'color', 'w');
    xlimtarg=[0.5:4:120.5];
    xticks(xlimtarg)
    xticklabels(bins(1:4:end))
    xlim([16.5 60.5])
    ylim([0.5 size(LickDayMat_Targets{ii}, 1)+0.5])
    title('Foil trials')
    
    hold all
    h3=subplot(1,10,7);
    h3.Position=[0.63 0.11 0.05 0.81];
    hold all
    plot(HitsR{ii},[1:1:size(LickDayMat_Foils{ii}, 1)], 'ko-', 'linewidth', 2) % reinforced hits
    plot(HitsP{ii},[1:1:size(LickDayMat_Foils{ii}, 1)], 'o-', 'color', [141, 141, 140]./255, 'linewidth', 2) %probe hits
    xlabel(['Hit proportion'])
    xlim([0 1])
    ylim([0.5 size(LickDayMat_Targets{ii}, 1)+0.5])
    title('Hits')
    
    h4=subplot(1,10,8);
    h4.Position=[0.72 0.11 0.05 0.81];
    hold all
    plot(FalseAlarmsR{ii},[1:1:size(LickDayMat_Foils{ii}, 1)], 'ko-', 'linewidth', 2) %reinforced fa
    plot(FalseAlarmsP{ii},[1:1:size(LickDayMat_Foils{ii}, 1)], 'o-', 'color', [141, 141, 140]./255, 'linewidth', 2) %probe fa
    xlabel(['False alarm proportion'])
    xlim([0 1])
    ylim([0.5 size(LickDayMat_Targets{ii}, 1)+0.5])
    title('False alarms')
    
    h5=subplot(1,10,10);
    h5.Position=[0.81 0.11 0.05 0.81];
    hold all
    plot(dprimeR{ii},[1:1:size(LickDayMat_Foils{ii}, 1)], 'ko-', 'linewidth', 2) %dprime reinforced
    plot(dprimePtouse{ii},[1:1:size(LickDayMat_Foils{ii}, 1)], 'o-', 'color', [141, 141, 140]./255, 'linewidth', 2) %dprime probe
    ylim([0.5 size(LickDayMat_Targets{ii}, 1)+0.5])
    xlabel(['d prime'])
    title('d-prime')
end

%% Mean licks per day per group
cooolor=[33, 158, 188; 255, 183, 3]./255;

for di=1:12
    wt_lickT=[];
    wt_lickF=[];
    mut_lickT=[];
    mut_lickF=[];
    for ii=1:length(data)
        
        if ismember(ii, wt)
            wt_lickT=[wt_lickT; LickDayMat_Targets{ii}(di,:)];
            wt_lickF=[wt_lickF; LickDayMat_Foils{ii}(di,:)];
%         elseif ismember(ii, mut)
%             mut_lickT=[mut_lickT; LickDayMat_Targets{ii}(di,:)];
%             mut_lickF=[mut_lickF; LickDayMat_Foils{ii}(di,:)];
        end
    end
    
    
    figure(1023)
    subplot(6,2,di)
    hold all
    shadedErrorBar(1:120,nanmean(wt_lickT),nanstd(wt_lickT)/sqrt(length(wt)), 'lineprops', {'-','color', cooolor(1,:), 'linewidth', 2})
%     shadedErrorBar(1:120,nanmean(mut_lickT),nanstd(mut_lickT)/sqrt(length(mut)), 'lineprops', {'-','color', cooolor(2,:), 'linewidth', 2})
    %plot(nanmean(wt_lickT), 'b')
    %plot(nanmean(mut_lickT), 'k')
    title(['Target licks day ' num2str(di)])
    ylim([0 80])
    
    figure(1024)
    subplot(6,2,di)
    hold all
    shadedErrorBar(1:120,nanmean(wt_lickF),nanstd(wt_lickF)/sqrt(length(wt)), 'lineprops', {'-','color', cooolor(1,:), 'linewidth', 2})
%     shadedErrorBar(1:120,nanmean(mut_lickF),nanstd(mut_lickF)/sqrt(length(mut)), 'lineprops', {'-','color', cooolor(2,:), 'linewidth', 2})
    %plot(nanmean(wt_lickT), 'b')
    %plot(nanmean(mut_lickT), 'k')
    title(['Foil licks day ' num2str(di)])
    ylim([0 40])
end