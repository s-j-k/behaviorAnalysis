%For figure making after pre-proccessed data
%S.M. 12/03/2019
%S.M. revised 04/40/2019
clear all
addpath('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Behavior\DSP4\')
addpath('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Matlab_functions\Outside');

datafrw=load('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Behavior\DSP4\Analysis\data_GNG_Normal_20210415_085535.mat');
datarev=load('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Behavior\DSP4\Analysis\data_GNG_Reverse_20210415_085319.mat');
data=datarev.data;

%% Find mice to use
for ii=1:length(data)
    dob(ii,:)=data(ii).dob;
    clear eini eend
    eini=num2str(data(ii).expday{1,1});
    eend=num2str(data(ii).expday{1,end});
    expini(ii,:)=[eini(5:6) '/' eini(7:8) '/' eini(1:4)];
    expend(ii,:)=[eend(5:6) '/' eend(7:8) '/' eend(1:4)];
    iniage(ii) = daysact(dob(ii,:), expini(ii,:))/30;
    endage(ii) = daysact(dob(ii,:), expend(ii,:))/30;
    tonego(ii,:)=unique(cell2mat(data(ii).tonego));
    tonenogo(ii,:)=unique(cell2mat(data(ii).tonenogo));
end 
meanage=(iniage+endage)/2;

%% To plot rasterplots of licks per animal per day
clear Licks LickMat Probetrials Reinforcedtrials BigMat

for ii=1:length(data)
    mouse{ii}=data(ii).mouse;
    gen(ii)=data(ii).genotype;
    trdays(ii)=length(data(ii).expday);
    sex(ii)=data(ii).sex;
    Licks=data(ii).lickstime; %data(ii).lickstime; lickstimefirst
    BigMat{ii}(1,:)=horzcat(data(ii).trialhist{:});%creates matrix with all info per trials not days for type of trial
    BigMat{ii}(2,:)=horzcat(data(ii).trialresp{:});%creates matrix with all info per trials not days for response
    dp_lucid{ii}=NaN(length(data(ii).expday),500);
    lucid_H{ii}=NaN(length(data(ii).expday),500);
    lucid_FA{ii}=NaN(length(data(ii).expday),500);
    
    %to calculate Probe of only first 10 trials
    trialsProbe=find(BigMat{ii}(1,:)>2);
    trialsProberesp=BigMat{ii}(2,trialsProbe);
    ini_winip=[1:20:1000];
    end_winip=[10:20:1000];
    for ppi=1:length(trialsProbe)/20
        probeblock=[ini_winip(ppi):end_winip(ppi)];
        ppiH=sum(trialsProberesp(probeblock)==1);
        ppiM=sum(trialsProberesp(probeblock)==2);
        ppiCR=sum(trialsProberesp(probeblock)==3);
        ppiFA=sum(trialsProberesp(probeblock)==4);   
        [dp_P10{ii}(ppi), c_P10{ii}(ppi)]=dprime_criterion_calc(ppiH, ppiM, ppiCR, ppiFA);
        ppi10_dp_P{ii}(ppi)=max(dp_P10{ii}(1:ppi));
        ppi10_xtrialP{ii}(ppi)=nanmean(trialsProbe(probeblock));
    end
    
    for di=1:length(data(ii).expday)
        ProbeTrials=find(data(ii).trialhist{1,di}>=3);
        ProbeTTargets=find(data(ii).trialhist{1,di}==3);
        ProbeTFoils=find(data(ii).trialhist{1,di}==4);
        ReinforcedTrials=find(data(ii).trialhist{1,di}<=2);
        Licksday=Licks{1,di};
        Targets=find(data(ii).trialhist{1,di}==1);
        Foils=find(data(ii).trialhist{1,di}==2);
        tlicksR{ii}(di,:)=data(ii).licksnum{1,di}(1);
        tlicksP{ii}(di,:)=data(ii).licksnum{1,di}(3);
        flicksR{ii}(di,:)=data(ii).licksnum{1,di}(2);
        flicksP{ii}(di,:)=data(ii).licksnum{1,di}(4);
        clear TrialResponse TtrialHist
        TrialResponse=data(ii).trialresp{1,di};
        TtrialHist=data(ii).trialhist{1,di};
        firstlickonset_RT{ii}(di,:)=nanmean(data(ii).firstlickonset{1,di}.RT);
        firstlickonset_RF{ii}(di,:)=nanmean(data(ii).firstlickonset{1,di}.RF);
        
        confidenceint=[0.2 0.5];
        RTconfident{ii}(di,:)=(length(find(data(ii).firstlickonset{1,di}.RT>confidenceint(1) & data(ii).firstlickonset{1,di}.RT<confidenceint(2)))/length(data(ii).firstlickonset{1,di}.RT));
        RFconfident{ii}(di,:)=(length(find(data(ii).firstlickonset{1,di}.RF>confidenceint(1) & data(ii).firstlickonset{1,di}.RF<confidenceint(2)))/length(data(ii).firstlickonset{1,di}.RF));

%         %%
%         %block analysis
% %         if length(ReinforcedTrials)>=240
% %             block1=[1:50];
% %             block2=[91:140];
% %             block3=[181:230];
% %         else
%             block1=[1:50];
%             block2=[91:140];
%             block3=[181:230];
% %         end
%         
%         if length(TtrialHist)<max(block3)
%             dp_rblock{ii}(di,1:3)=NaN;
%             c_rblock{ii}(di,1:3)=NaN;
%             lucidd{ii}(di,:)=NaN;
%         else
%             Rblock1=ReinforcedTrials(block1);
%             Rblock2=ReinforcedTrials(block2);
%             Rblock3=ReinforcedTrials(block3);
%             blocks=[Rblock1; Rblock2; Rblock3];
%             
%             for bli=1:size(blocks,1)
%                 Rhitsb=sum(TrialResponse(blocks(bli,:))==1);
%                 Rmissb=sum(TrialResponse(blocks(bli,:))==2);
%                 Rcrb=sum(TrialResponse(blocks(bli,:))==3);
%                 Rfab=sum(TrialResponse(blocks(bli,:))==4);
%                 [dp_rblock{ii}(di,bli), c_rblock{ii}(di,bli)]=dprime_criterion_calc(Rhitsb, Rmissb, Rcrb, Rfab);
%             end
%             
%             %for lucidity periods
            slidwindow=[1:30];
            maxdprime=data(ii).dprime{1,di}(1);
            winincrease=1;
% 
            for lui=1:(length(ReinforcedTrials)-length(slidwindow))/winincrease
                Lhitsb=sum(TrialResponse(ReinforcedTrials(slidwindow))==1);
                Lmissb=sum(TrialResponse(ReinforcedTrials(slidwindow))==2);
                Lcrb=sum(TrialResponse(ReinforcedTrials(slidwindow))==3);
                Lfab=sum(TrialResponse(ReinforcedTrials(slidwindow))==4);
                [dp_lucid{ii}(di,lui), c_lucid{ii}(di,lui)]=dprime_criterion_calc(Lhitsb, Lmissb, Lcrb, Lfab);
                lucid_H{ii}(di,lui)=Lhitsb/(Lhitsb+Lmissb);
                lucid_FA{ii}(di,lui)=Lfab/(Lfab+Lcrb);
                slidwindow=slidwindow+winincrease;
            end
% %             maxd=nanmean(dp_lucid{ii}(di,1:lui));
% %             figure(ii*winincrease);  
% %             set(gcf, 'Position', get(0, 'Screensize'))
% %             subplot(6,6,di)
% %             hold all
% %             plot(dp_lucid{ii}(di,:), 'k-', 'linewidth', 2)
% %             line([0 300], [maxdprime maxdprime], 'color', 'r')
% %             line([0 300], [maxd maxd], 'color', 'g')
% %             line([0 300], [dprimePtoplot(di,ii) dprimePtoplot(di,ii)], 'color', 'b')
% %             line([0 300], [max(dprimePtoplot(:,ii)) max(dprimePtoplot(:,ii))], 'color', 'm')
% %             ylabel('d''')
% %             xlabel(['Session window (' num2str(length(slidwindow)) ' trials)'])
% %             title([mouse{ii} filesep num2str(di)])
% %             ylim([-2 5])
%         end
%%        
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
        
        
%         if di==1
%             dprimeRtouse{ii}(di,:)=dprimeR{ii}(di,:);
%         elseif di>=2
%             if dprimeR{ii}(di,:)>=dprimeRtouse{ii}(di-1,:)
%                 dprimeRtouse{ii}(di,:)=dprimeR{ii}(di,:);
%             else
%                 dprimeRtouse{ii}(di,:)=max(dprimeRtouse{ii});
%             end
%         end
        
        %to use only higher dprimes during probe 
        if isnan(dprimeP{ii}(di,:))
            dprimePtouse{ii}(di,:)=NaN;
        elseif ~isnan(dprimeP{ii}(di,:))
            if sum(isnan(dprimePtouse{ii}))==length(dprimePtouse{ii}) %if this is the first numeric d prime to be added
                dprimePtouse{ii}(di,:)=dprimeP{ii}(di,:);
            elseif dprimeP{ii}(di,:)>=max(dprimePtouse{ii}) %if other max dprimes exist
                dprimePtouse{ii}(di,:)=dprimeP{ii}(di,:);
            else
                dprimePtouse{ii}(di,:)=max(dprimePtouse{ii});
            end
        end
        
    end
end

%% FORMAT matrices for plotting
dprimeRtoplot=NaN(40, length(data));
dprimePtoplot=NaN(40, length(data));
hitRtoplot=NaN(40, length(data));
faRtoplot=NaN(40, length(data));
hitPtoplot=NaN(40, length(data));
faPtoplot=NaN(40, length(data));
tlicksRtoplot=NaN(40, length(data));
tlicksPtoplot=NaN(40, length(data));
flicksRtoplot=NaN(40, length(data));
flicksPtoplot=NaN(40, length(data));
RTlicktimetoplot=NaN(40, length(data));
RFlicktimetoplot=NaN(40, length(data));
RTconfidenttoplot=NaN(40, length(data));
RFconfidenttoplot=NaN(40, length(data));

% limit dp probe performance

mindprobedp=0.75;

for ii=1:length(data)
dprimeRtoplot(1:length(dprimeR{ii}),ii)=dprimeR{ii};
dprimePtoplot(1:length(dprimePtouse{ii}),ii)=dprimePtouse{ii};
hitRtoplot(1:length(HitsR{ii}),ii)=HitsR{ii};
hitPtoplot(1:length(HitsP{ii}),ii)=HitsP{ii};
faRtoplot(1:length(FalseAlarmsR{ii}),ii)=FalseAlarmsR{ii};
faPtoplot(1:length(FalseAlarmsP{ii}),ii)=FalseAlarmsP{ii};
tlicksRtoplot(1:length(tlicksR{ii}), ii)=tlicksR{ii};
tlicksPtoplot(1:length(tlicksP{ii}), ii)=tlicksP{ii};
flicksRtoplot(1:length(flicksR{ii}), ii)=flicksR{ii};
flicksPtoplot(1:length(flicksP{ii}), ii)=flicksP{ii};
RTlicktimetoplot(1:length(firstlickonset_RT{ii}), ii)=firstlickonset_RT{ii};
RFlicktimetoplot(1:length(firstlickonset_RF{ii}), ii)=firstlickonset_RF{ii};
RTconfidenttoplot(1:length(RTconfident{ii}), ii)=RTconfident{ii};
RFconfidenttoplot(1:length(RFconfident{ii}), ii)=RFconfident{ii};


if sum(find(dprimeRtoplot(:,ii)>=mindprobedp))>=1
    learnedR(ii)=1;
else
    learnedR(ii)=NaN;
end 

if sum(find(dprimePtoplot(:,ii)>=mindprobedp))>=1
    learnedP(ii)=1;
else
    learnedP(ii)=NaN;
end
end
 
%% FIND MICE OF INTEREST
gen=[data(1:28).genotype];
% ALL mice
WTtouseyoung=find(gen==1);
% MUTtouseyoung=find(gen==2); %[3 4 7 8 10 12]; %
MUTtouseyoung=find(gen==2); %[3 4 7 8 10 12]; %
WTtouseyoungN=WTtouseyoung(1:2);
WTtouseyoungR=WTtouseyoung(3:4);
MUTtouseyoungN=MUTtouseyoung(1:3);
MUTtouseyoungR=MUTtouseyoung(4:6);

% WTtouse=intersect(WTtouse, find(learnedP==1));
% MUTtouse=intersect(MUTtouse, find(learnedP==1));
% 
% WT_mice=mouse(WTtouse); %for wt
% MUT_mice=mouse(MUTtouse); %for wt
% 
% %older vs younger mice
% WTtouseyoung=intersect(WTtouse, find(meanage<9));
% WTtouseold=intersect(WTtouse, find(meanage>9));
% MUTtouseyoung=intersect(MUTtouse, find(meanage<9));
% MUTtouseold=intersect(MUTtouse, find(meanage>9));
% 
% AllWTold=intersect(find(gen==1), find(meanage>9));
% AllMUTold=intersect(find(gen==2), find(meanage>9));

%% Sliding window d' every 100 trials for R and every 20 for P
clear dp_slideR pHit pFA
colorsmat2=[96, 211, 148; 238, 96, 85]./255; %[[155, 204, 85]./255; [228, 87, 46]./255; [179, 201, 174]./255; [239, 149, 137]./255];
dp_slideR=NaN(length(data),8000);
c_slideR=NaN(length(data),8000);
pHit=NaN(length(data),200);
pFA=NaN(length(data),200);

winlength=30;
winincrease=2;
winlength2=100;
winincrease2=100;

for ii=1:length(data)
    clear ReinforcedTrials TtrialHist TrialResponse
    ReinforcedTrials=find(BigMat{ii}(1,:)<=2);
    TtrialHist=BigMat{ii}(1,ReinforcedTrials);
    TrialResponse=BigMat{ii}(2,ReinforcedTrials);
    slidwindow=[1:winlength];
    for sli=1:(length(TtrialHist)-length(slidwindow))/winincrease
        Lhitsb=sum(TrialResponse(slidwindow)==1);
        Lmissb=sum(TrialResponse(slidwindow)==2);
        Lcrb=sum(TrialResponse(slidwindow)==3);
        Lfab=sum(TrialResponse(slidwindow)==4);
        [dp_slideR(ii,sli), c_slideR(ii,sli)]=dprime_criterion_calc(Lhitsb, Lmissb, Lcrb, Lfab);
        slidwindow=slidwindow+winincrease;
    end
    
    slidwindow2=[1:winlength2];
    for wi=1:round(length(TtrialHist)/winincrease2)-1
       pHit(ii,wi)=sum(TrialResponse(slidwindow2)==1)/sum(TtrialHist(slidwindow2)==1);
       pFA(ii,wi)=sum(TrialResponse(slidwindow2)==4)/sum(TtrialHist(slidwindow2)==2);
       slidwindow2=slidwindow2+winincrease2;
    end 
end

figure
i=0;
for ii=MUTtouseyoung
    i=i+1;
    subplot(3,4,i)
    hold all
    plot(dp_slideR(ii,:)', 'k-')
%     plot(nanmean(dp_slideR(WTtouseyoung,:))', 'k-', 'linewidth', 2)
%     title([char(mouse(ii)) filesep num2str(winlength) filesep num2str(winincrease)])
    ylabel('d''')
    ylim([-2 5])
    xlabel('Bin trial #')
    line([0 5410], [max(dprimePtoplot(:,ii)) max(dprimePtoplot(:,ii))], 'color', 'r')
end

figure
subplot(1,2,1)
hold all
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
% for ii=WTtouseold
    plot([100:100:8000], nanmean(pHit(WTtouseyoung,1:80))', '-', 'color', colorsmat2(1,:), 'linewidth', 2); %, 'Markersize', 15)
    plot([100:100:8000], nanmean(pFA(WTtouseyoung,1:80))', '-', 'color', colorsmat2(2,:), 'linewidth', 2); %, 'Markersize', 15)
    ylim([0 1]); xlim([1 8000]);
    xlabel('Trial #')
    ylabel('Action rate')
    title('Reinforced WR')
% end
subplot(1,2,2)
hold all
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
% for ii=MUTtouseold
    plot([100:100:8000], nanmean(pHit(MUTtouseyoung,1:80))', '-', 'color', colorsmat2(1,:), 'linewidth', 2); %, 'Markersize', 15)
    plot([100:100:8000], nanmean(pFA(MUTtouseyoung,1:80))', '-', 'color', colorsmat2(2,:), 'linewidth', 2); %, 'Markersize', 15)
    ylim([0 1]); xlim([1 8000]);
    xlabel('Trial #')
    ylabel('Action rate')
    title('Reinforced CA')
 
    


%% D PRIME ANALYSIS PER TRIALS GOOD
colors=[[187, 186, 183]; [255, 195, 185]]./255;

winlength3=100;
winincrease3=100;
winlength4=20;
winincrease4=20;
maxsizemat=500;

clear dp_sli400R dp_sli400P dp_slitouseP
xtrialsP=NaN(length(data),maxsizemat);
xtrialsR=NaN(length(data),maxsizemat);
dp_slitouseP=NaN(length(data),maxsizemat);
dp_sli100P=NaN(length(data),maxsizemat);
dp_sli100R=NaN(length(data),maxsizemat);
c_sli100R=NaN(length(data),maxsizemat);
pHsli100R=NaN(length(data),maxsizemat);
pFAsli100R=NaN(length(data),maxsizemat);
pHsli100P=NaN(length(data),maxsizemat);
pFAsli100P=NaN(length(data),maxsizemat);
maxdpsli100R=NaN(length(data),maxsizemat);
maxcsli100R=NaN(length(data),maxsizemat);
pH_sli100R=NaN(length(data),maxsizemat);
pFA_sli100R=NaN(length(data),maxsizemat);

for ii=1:length(data)
    clear ReinforcedTrials ProbeTrials TtrialHistR TtrialHistP TrialResponseR TrialResponseP
    ReinforcedTrials=find(BigMat{ii}(1,:)<=2);
    ProbeTrials=find(BigMat{ii}(1,:)>2);
    TtrialHistR=BigMat{ii}(1,ReinforcedTrials);
    TrialResponseR=BigMat{ii}(2,ReinforcedTrials);
    TtrialHistP=BigMat{ii}(1,ProbeTrials);
    TrialResponseP=BigMat{ii}(2,ProbeTrials);
    NumTones(ii,:)=[sum(TtrialHistR==1)/length(TtrialHistR) sum(TtrialHistR==2)/length(TtrialHistR)];
    sliwi=0; %for reinforced
    for sli=1:((max(ReinforcedTrials)-winincrease3)/winincrease3)
        clear til
        til=find(ReinforcedTrials>sliwi & ReinforcedTrials<=sliwi+winincrease3);
        tillen{ii}(sli,:)=length(til);
        xtrialsR(ii,sli)=round(nanmean(ReinforcedTrials(til)));
        [maxdpsli100R(ii,sli) maxcsli100R(ii,sli)]=dprime_criterion_calc(sum(TtrialHistR(til)==1), 0, sum(TtrialHistR(til)==2), 0);
        Rhitssli=sum(TrialResponseR(til)==1);
        Rmisssli=sum(TrialResponseR(til)==2);
        Rcrsli=sum(TrialResponseR(til)==3);
        Rfasli=sum(TrialResponseR(til)==4);
        [dp_sli100R(ii,sli), c_sli100R(ii,sli)]=dprime_criterion_calc(Rhitssli, Rmisssli, Rcrsli, Rfasli);
        [bla1,bla2,pH_sli100R(ii,sli), pFA_sli100R(ii,sli)]=dprime_criterion_pH_pFA_calc(Rhitssli, Rmisssli, Rcrsli, Rfasli);
        pHsli100R(ii,sli)=Rhitssli/(Rhitssli+Rmisssli);
        pFAsli100R(ii,sli)=Rfasli/(Rfasli+Rcrsli);
        sliwi=sliwi+winincrease3;
    end
    rrr(ii,:)=min(tillen{ii});
     
    sliwi=0;
    slidwindow4=[1:winlength4]; %for probe
    for sli=1:length(TtrialHistP)/winincrease4
        clear til        
        xtrialsP(ii,sli)=round(nanmean(ProbeTrials(slidwindow4)));
        Phitssli=sum(TrialResponseP(slidwindow4)==1);
        Pmisssli=sum(TrialResponseP(slidwindow4)==2);
        Pcrsli=sum(TrialResponseP(slidwindow4)==3);
        Pfasli=sum(TrialResponseP(slidwindow4)==4);
        [dp_sli100P(ii,sli), c_sli100P(ii,sli)]=dprime_criterion_calc(Phitssli, Pmisssli, Pcrsli, Pfasli);
        pHsli100P(ii,sli)=Phitssli/(Phitssli+Pmisssli);
        pFAsli100P(ii,sli)=Pfasli/(Pfasli+Pcrsli);
        slidwindow4=slidwindow4+winincrease4;
        if sli>1
            if dp_sli100P(ii,sli)>=max(dp_sli100P(ii,:))
                dp_slitouseP(ii,sli)=dp_sli100P(ii,sli);
            else
                dp_slitouseP(ii,sli)=max(dp_sli100P(ii,:));
            end
        else
            dp_slitouseP(ii,sli)=dp_sli100P(ii,sli);
        end
    end
    
    %lucidity sliding window
%     slidwindow=[1:20];
%     winincrease=1;
%     for lui=1:(length(ReinforcedTrials)-length(slidwindow))/winincrease
%         Lhitsb=sum(TrialResponse(ReinforcedTrials(slidwindow))==1);
%         Lmissb=sum(TrialResponse(ReinforcedTrials(slidwindow))==2);
%         Lcrb=sum(TrialResponse(ReinforcedTrials(slidwindow))==3);
%         Lfab=sum(TrialResponse(ReinforcedTrials(slidwindow))==4);
%         [dp_lucid_sli100R{ii}(di,lui), c_lucid_sli100R{ii}(di,lui)]=dprime_criterion_calc(Lhitsb, Lmissb, Lcrb, Lfab);
%         lucid_sli100R_H{ii}(di,lui)=Lhitsb/(Lhitsb+Lmissb);
%         lucid_sli100R_FA{ii}(di,lui)=Lfab/(Lfab+Lcrb);
%         slidwindow=slidwindow+winincrease;
%     end
end

%% normalized max dprime
plotscolor=[124, 77, 234]./255;
wtmice=WTtouseyoung;
mutmice=MUTtouseyoung;
maxday=500;

figure
hold all
for ii=wtmice
plot([1:maxday]*100,(dp_sli100R(ii,1:maxday)*100)./maxdpsli100R(ii,1:maxday), 'k-', 'linewidth', 2)
ylabel('Normalized performance')
xlabel('Trials')
text([200],[80], 'saline', 'Color', 'k','FontSize', 12, 'fontname','arial')
text([200],[70], 'dsp4', 'Color', plotscolor,'FontSize', 12, 'fontname','arial')
end 
for ii=mutmice
plot([1:maxday]*100,(dp_sli100R(ii,1:maxday)*100)./maxdpsli100R(ii,1:maxday), '-', 'color', plotscolor, 'linewidth', 2)
end 

for ii=1:length(mouse)
    maxdpsli100Rnorm(ii,:)=max((dp_sli100R(ii,1:maxday)*100)./maxdpsli100R(ii,1:maxday));
    maxdpsli100Pnorm(ii,:)=max((dp_slitouseP(ii,1:maxday)*100)./3.2897);
end

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, maxdpsli100Rnorm(wtmice)', [0,0,0],6,1);
paper_dotbarplot_modif(1, maxdpsli100Rnorm(mutmice)', plotscolor,6,1);
paper_dotbarplot_modif(1.8, maxdpsli100Pnorm(wtmice)', [0,0,0],6,1);
paper_dotbarplot_modif(2.3, maxdpsli100Pnorm(mutmice)', plotscolor,6,1);
ylabel('Normalized performance')
set(gca, 'Xtick', [0.5 1 1.8 2.3], 'Xticklabel', {'saline', 'dsp4'})
ylim([0 100])
line([0 2.8],[80 80], 'color', 'k', 'linestyle', ':')
xlim([0 2.8])
title('Reinforced                             Probe')

%% ROC curves
wtmice=WTtouseyoung;
mutmice=MUTtouseyoung;



figure
hold all
for ii=1%wtmice
%     ROCx=
    plot(pFA_sli100R(ii,:), pH_sli100R(ii,:), 'bo')
    ylabel('Hit probability')
    xlabel('FA probability')
    xlim([0 1]); ylim([0 1]);
    grid on
    line([0 1],[0 1],'color', 'k', 'linestyle', ':')
end
for ii=mutmice
    plot(pFA_sli100R(ii,:), pH_sli100R(ii,:), 'ko')
    xlim([0 1]); ylim([0 1]);
    grid on
    line([0 1],[0 1],'color', 'k', 'linestyle', ':')
end

%% ROC curves
wtmice=WTtouseyoung;
mutmice=MUTtouseyoung;


clear ROCx_wt ROCy_wt ROCx_mut ROCy_mut

i=0;
for ii=wtmice
    i=i+1;
    ROCx_wt(i,:)=pFA_sli100R(ii,:);
    ROCy_wt(i,:)=pH_sli100R(ii,:);
end

i=0;
for ii=mutmice
    i=i+1;
    ROCx_mut(i,:)=pFA_sli100R(ii,:);
    ROCy_mut(i,:)=pH_sli100R(ii,:); 
end


figure
hold all
plot(c_sli100R(wtmice,:), ROCy_wt, 'go')
plot(c_sli100R(wtmice,:), ROCx_wt, 'ro')
xlabel('Criterion')
ylabel('Proportion')
title('saline')
xlim([-2.5 2.5])


figure
hold all
plot(c_sli100R(mutmice,:), ROCy_mut, 'go')
plot(c_sli100R(mutmice,:), ROCx_mut, 'ro')
xlabel('Criterion')
ylabel('Proportion')
title('dsp4')
xlim([-2.5 2.5])

figure
hold all
plot(nanmean(ROCx_mut), nanmean(ROCy_mut), 'bo')
plot(nanmean(ROCx_wt), nanmean(ROCy_wt), 'ko')
xlabel('FA probability')
ylabel('Hit probability')
grid on
line([0 1],[0 1],'color', 'k', 'linestyle', ':')
axis square

%% Dprime vs Criterion comparison
wtmice=WTtouseyoung;
mutmice=MUTtouseyoung;
usedwinlength=100;
maxtrialstoplot=4000;

figure
hold all
i=0;
for ii=wtmice
    i=i+1;
    subplot(4,4,i)
    hold all
    yyaxis left
    plot([1:usedwinlength:maxtrialstoplot], dp_sli100R(ii,1:maxtrialstoplot/usedwinlength), 'b')
    ylabel('d'''); xlim([1 maxtrialstoplot]); ylim([-2 5]); 
    yyaxis right
    plot([1:usedwinlength:maxtrialstoplot], -1*(c_sli100R(ii,1:maxtrialstoplot/usedwinlength)), 'r')
    line([1 maxtrialstoplot],[-1.2 -1.2], 'color', 'k', 'linestyle', ':')
    line([1 maxtrialstoplot],[1.2 1.2], 'color', 'k', 'linestyle', ':')
    ylabel('-c'); xlim([1 maxtrialstoplot]); ylim([-3 3]);
    title(['\color{Blue}' char(mouse(ii))])
end
i=8;
for ii=mutmice
    i=i+1;
    subplot(4,4,i)
    hold all
    yyaxis left
    plot([1:usedwinlength:maxtrialstoplot], dp_sli100R(ii,1:maxtrialstoplot/usedwinlength), 'b')
    ylabel('d'''); xlim([1 maxtrialstoplot]); ylim([-2 5]);
    yyaxis right
    plot([1:usedwinlength:maxtrialstoplot],-1*(c_sli100R(ii,1:maxtrialstoplot/usedwinlength)), 'r')
    line([1 maxtrialstoplot],[-1.2 -1.2], 'color', 'k', 'linestyle', ':')
    line([1 maxtrialstoplot],[1.2 1.2], 'color', 'k', 'linestyle', ':')
    ylabel('-c'); xlim([1 maxtrialstoplot]); ylim([-3 3]);
    title([char(mouse(ii))])
end


%% MAX ever dp in reinforced and probe context
colors=[0, 184, 245; 254, 221, 11]./255;
wtmice=WTtouseyoung; %AllWTold; %WTtouseold;
mutmice=MUTtouseyoung;%AllMUTold; %MUTtouseold;

%per trial
for ii=1:length(data)
    maxdpR(ii, :)=max(dp_sli100R(ii,:));
    maxdpP(ii, :)=max(dp_slitouseP(ii,:));
end 

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, maxdpR(wtmice)', colors(1,:),6,1);
paper_dotbarplot_modif(1, maxdpR(mutmice)', colors(2,:),6,1);
paper_dotbarplot_modif(1.8, maxdpP(wtmice)', colors(1,:),6,1);
paper_dotbarplot_modif(2.3, maxdpP(mutmice)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1 1.8 2.3], 'Xticklabel', {'WR', 'CA'})
ylabel('Max d'' overall')
% ylim([0 1.05])
xlim([0 2.8])
pval=mynormandttest(maxdpR(wtmice), maxdpR(mutmice))
pval=mynormandttest(maxdpP(wtmice), maxdpP(mutmice))

%per day
for ii=1:length(data)
    maxdpR_day(ii, :)=max(dprimeR{ii});
    maxdpP_day(ii, :)=max(dprimeP{ii});
end 

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, maxdpR_day(wtmice)', colors(1,:),6,1);
paper_dotbarplot_modif(1, maxdpR_day(mutmice)', colors(2,:),6,1);
paper_dotbarplot_modif(1.8, maxdpP_day(wtmice)', colors(1,:),6,1);
paper_dotbarplot_modif(2.3, maxdpP_day(mutmice)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1 1.8 2.3], 'Xticklabel', {'WR', 'CA'})
ylabel('Max d'' overall per day')
% ylim([0 1.05])
xlim([0 2.8])
pval=mynormandttest(maxdpR_day(wtmice), maxdpR_day(mutmice))
pval=mynormandttest(maxdpP_day(wtmice), maxdpP_day(mutmice))
%% 
%% FIGURE 2: mean PH pFA action rates for Reinforced and Probe
wtmice=[WTtouseyoungN] ;% WTtouseold]; %WTtouseyoung; %WTtouseold;
mutmice=[MUTtouseyoungN]; % MUTtouseold]; %MUTtouseyoung; %MUTtouseold;

addpath('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Matlab_functions\Plotting\functions_plots\')
clear xp yp semp
maxdtop=40; 

colorsmat2=[96, 211, 148; 242, 143, 59]./255; % 238, 96, 85 . [[155, 204, 85]./255; [228, 87, 46]./255; [179, 201, 174]./255; [239, 149, 137]./255];

yp=[nanmean(pHsli100R(wtmice,1:maxdtop));nanmean(pFAsli100R(wtmice,1:maxdtop));nanmean(pHsli100R(mutmice,1:maxdtop));nanmean(pFAsli100R(mutmice,1:maxdtop));...
    nanmean(pHsli100P(wtmice,1:maxdtop));nanmean(pFAsli100P(wtmice,1:maxdtop));nanmean(pHsli100P(mutmice,1:maxdtop));nanmean(pFAsli100P(mutmice,1:maxdtop))];
xp=[round(nanmean(xtrialsR(wtmice, 1:maxdtop)));round(nanmean(xtrialsR(wtmice, 1:maxdtop)));round(nanmean(xtrialsR(mutmice, 1:maxdtop)));round(nanmean(xtrialsR(mutmice, 1:maxdtop)));...
    round(nanmean(xtrialsP(wtmice, 1:maxdtop)));round(nanmean(xtrialsP(wtmice, 1:maxdtop)));round(nanmean(xtrialsP(mutmice, 1:maxdtop)));round(nanmean(xtrialsP(mutmice, 1:maxdtop)))];
colp=[colorsmat2(1,:);colorsmat2(2,:);colorsmat2(1,:);colorsmat2(2,:);colorsmat2(1,:);colorsmat2(2,:);colorsmat2(1,:);colorsmat2(2,:)];

semp=[semcalc(pHsli100R(wtmice,1:maxdtop)); semcalc(pFAsli100R(wtmice,1:maxdtop)); semcalc(pHsli100R(mutmice,1:maxdtop)); semcalc(pFAsli100R(mutmice,1:maxdtop));...
      semcalc(pHsli100P(wtmice,1:maxdtop)); semcalc(pFAsli100P(wtmice,1:maxdtop)); semcalc(pHsli100P(mutmice,1:maxdtop)); semcalc(pFAsli100P(mutmice,1:maxdtop))]; 

mp=[mouse(wtmice) mouse(mutmice)];

titl=[{'Reinforced s','Reinforced dsp4','Probe s','Probe dsp4'}];  

figure
fi=1;
for si=1:4
    subplot(2,2,si); hold all
    set(gca,'TickDir','out','FontSize',12,'fontname','arial');
    title(char(titl(si))); ylabel('Action rate'); xlabel('Trials'); ylim([0 1.02]); xlim([0 maxdtop*100]);
    if si==1 || si==2 %for Reinforced
        shadedErrorBar(xp(fi,:),yp(fi,:),semp(fi,:), 'lineprops', {'-','color', colorsmat2(1,:), 'linewidth', 2})
        shadedErrorBar(xp(fi+1,:),yp(fi+1,:),semp(fi+1,:), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
    else %for Probe
        shadedErrorBar(xp(fi,1:12),yp(fi,1:12),semp(fi,1:12), 'lineprops', {'-','color', colorsmat2(1,:), 'linewidth', 2})
        shadedErrorBar(xp(fi+1,1:12),yp(fi+1,1:12),semp(fi+1,1:12), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
    end
    fi=fi+2;
    if si==1
        text([6100],[0.99], 'Target', 'Color', colorsmat2(1,:),'FontSize', 12, 'fontname','arial')
        text([6100], [0.1], 'Foil', 'Color', colorsmat2(2,:),'FontSize', 12, 'fontname','arial')
    end
end

%% FIGURE 3 NSF grant Kishore: mean PH pFA action rates for Reinforced and Probe
colorsmat2=[96, 211, 148; 242, 143, 59]./255; % 238, 96, 85 . [[155, 204, 85]./255; [228, 87, 46]./255; [179, 201, 174]./255; [239, 149, 137]./255];
addpath('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Matlab_functions\Plotting\functions_plots\')
maxdtop=30;
clear reinforcedRev probeRev

for pli=1:2
    
    clear wtmice mutmice
    if pli==1
        wtmice=[WTtouseyoungN] ;
        mutmice=[MUTtouseyoungN];
    elseif  pli==2
        wtmice=[WTtouseyoungR] ;
        mutmice=[MUTtouseyoungR];
    end
    
    clear xp yp semp
    yp=[nanmean(pHsli100R(wtmice,1:maxdtop));nanmean(pFAsli100R(wtmice,1:maxdtop));nanmean(pHsli100R(mutmice,1:maxdtop));nanmean(pFAsli100R(mutmice,1:maxdtop));...
        nanmean(pHsli100P(wtmice,1:maxdtop));nanmean(pFAsli100P(wtmice,1:maxdtop));nanmean(pHsli100P(mutmice,1:maxdtop));nanmean(pFAsli100P(mutmice,1:maxdtop))];
    xp=[round(nanmean(xtrialsR(wtmice, 1:maxdtop)));round(nanmean(xtrialsR(wtmice, 1:maxdtop)));round(nanmean(xtrialsR(mutmice, 1:maxdtop)));round(nanmean(xtrialsR(mutmice, 1:maxdtop)));...
        round(nanmean(xtrialsP(wtmice, 1:maxdtop)));round(nanmean(xtrialsP(wtmice, 1:maxdtop)));round(nanmean(xtrialsP(mutmice, 1:maxdtop)));round(nanmean(xtrialsP(mutmice, 1:maxdtop)))];
    colp=[colorsmat2(1,:);colorsmat2(2,:);colorsmat2(1,:);colorsmat2(2,:);colorsmat2(1,:);colorsmat2(2,:);colorsmat2(1,:);colorsmat2(2,:)];
    
    semp=[semcalc(pHsli100R(wtmice,1:maxdtop)); semcalc(pFAsli100R(wtmice,1:maxdtop)); semcalc(pHsli100R(mutmice,1:maxdtop)); semcalc(pFAsli100R(mutmice,1:maxdtop));...
        semcalc(pHsli100P(wtmice,1:maxdtop)); semcalc(pFAsli100P(wtmice,1:maxdtop)); semcalc(pHsli100P(mutmice,1:maxdtop)); semcalc(pFAsli100P(mutmice,1:maxdtop))];
    
    mp=[mouse(wtmice) mouse(mutmice)];
    
    titl=[{'Reinforced s','Reinforced dsp4','Probe s','Probe dsp4'}];
    
    figure(1032)
    fi=1:2;
    for si=1
%         subplot(1,2,si); hold all
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
        title(char(titl(si))); ylabel('Action rate'); xlabel('Trials'); ylim([0 1.02]); xlim([0 2*maxdtop*100]);
%         if si==1 || si==2 %for Reinforced
         if pli==1
%             shadedErrorBar(xp(fi,:),yp(fi,:),semp(fi,:), 'lineprops', {'-','color', colorsmat2(1,:), 'linewidth', 2})
%             shadedErrorBar(xp(fi+1,:),yp(fi+1,:),semp(fi+1,:), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
         elseif pli==2  
%             shadedErrorBar([2500]+xp(fi,:),yp(fi,:),semp(fi,:), 'lineprops', {'-','color', colorsmat2(1,:), 'linewidth', 2})
%             shadedErrorBar([2500]+xp(fi+1,:),yp(fi+1,:),semp(fi+1,:), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
%         else %for Probe
%             shadedErrorBar(xp(fi,1:12),yp(fi,1:12),semp(fi,1:12), 'lineprops', {'-','color', colorsmat2(1,:), 'linewidth', 2})
%             shadedErrorBar(xp(fi+1,1:12),yp(fi+1,1:12),semp(fi+1,1:12), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
        end
        fi=fi+2;
        if si==1
            text([5450],[0.99], 'Target', 'Color', colorsmat2(1,:),'FontSize', 12, 'fontname','arial')
            text([5450], [0.1], 'Foil', 'Color', colorsmat2(2,:),'FontSize', 12, 'fontname','arial')
        end
    end
    
    
% To store mean data for model
if pli==1
    reinforcedRev(1:25,1)=xp(1,1:25)'; %trials
    reinforcedRev(1:25,2)=yp(1,1:25)'; %hits
    reinforcedRev(1:25,3)=yp(2,1:25)'; %fa
    
    probeRev(1:7,1)=xp(5,1:7)'; %trials probe
    probeRev(1:7,2)=yp(5,1:7)'; %hits
    probeRev(1:7,3)=yp(6,1:7)'; %fa
    
elseif pli==2 
    reinforcedRev(26:55,1)=2500+xp(1,1:30)'; %trials
    reinforcedRev(26:55,2)=yp(1,1:30)'; %hits
    reinforcedRev(26:55,3)=yp(2,1:30)'; %fa
    
    probeRev(8:12,1)=2500+xp(5,1:5)'; %trials probe
    probeRev(8:12,2)=yp(5,1:5)'; %hits
    probeRev(8:12,3)=yp(6,1:5)'; %fa
end

end
line([2500 2500], [0 1], 'color', 'k', 'linestyle', ':', 'linewidth',1.2);
text([2550],[0.99], 'Reversal', 'Color', 'k','FontSize', 12, 'fontname','arial')
text([15],[0.1], 'n=2', 'Color', 'k','FontSize', 12, 'fontname','arial')






%% FIGURE 4SF grant Kishore: mean dprime for Reinforced 
colorsmat2=[96, 211, 148; 242, 143, 59]./255; % 238, 96, 85 . [[155, 204, 85]./255; [228, 87, 46]./255; [179, 201, 174]./255; [239, 149, 137]./255];
addpath('\\pbs-srv2.win.ad.jhu.edu\KishoreLab\Sharlen\Matlab_functions\Plotting\functions_plots\')
maxdtop=30;

for pli=1:2
    
    clear wtmice mutmice
    if pli==1
        wtmice=[WTtouseyoungN] ;
        mutmice=[MUTtouseyoungN];
    elseif  pli==2
        wtmice=[WTtouseyoungR] ;
        mutmice=[MUTtouseyoungR];
    end
    
    clear xp yp semp
    yp=[nanmean(dp_sli100R(wtmice,1:maxdtop));nanmean(dp_sli100R(mutmice,1:maxdtop))];
    xp=[round(nanmean(xtrialsR(wtmice, 1:maxdtop)));round(nanmean(xtrialsR(mutmice, 1:maxdtop)))];    
    semp=[semcalc(dp_sli100R(wtmice,1:maxdtop)); semcalc(dp_sli100R(mutmice,1:maxdtop))];
    mp=[mouse(wtmice) mouse(mutmice)];
        
    figure(1035)
    fi=1;
    for si=1
%         subplot(1,2,si); hold all
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
%         title(char(titl(si))); 
        ylabel('d prime'); xlabel('Trials'); %lim([0 1.02]); 
        xlim([0 2*maxdtop*100]);
%         if si==1 || si==2 %for Reinforced
         if pli==1
            shadedErrorBar(xp(fi,:),yp(fi,:),semp(fi,:), 'lineprops', {'-','color', 'k', 'linewidth', 2})
%             shadedErrorBar(xp(fi+1,:),yp(fi+1,:),semp(fi+1,:), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
         elseif pli==2  
            shadedErrorBar([2500]+xp(fi,:),yp(fi,:),semp(fi,:), 'lineprops', {'-','color', 'k', 'linewidth', 2})
%             shadedErrorBar([2500]+xp(fi+1,:),yp(fi+1,:),semp(fi+1,:), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
%         else %for Probe
%             shadedErrorBar(xp(fi,1:12),yp(fi,1:12),semp(fi,1:12), 'lineprops', {'-','color', colorsmat2(1,:), 'linewidth', 2})
%             shadedErrorBar(xp(fi+1,1:12),yp(fi+1,1:12),semp(fi+1,1:12), 'lineprops', {'-','color', colorsmat2(2,:), 'linewidth', 2})
        end
        fi=fi+2;
%         if si==1
%             text([5450],[0.99], 'Target', 'Color', colorsmat2(1,:),'FontSize', 12, 'fontname','arial')
%             text([5450], [0.1], 'Foil', 'Color', colorsmat2(2,:),'FontSize', 12, 'fontname','arial')
%         end
    end
end
line([2500 2500], [-2 5], 'color', 'k', 'linestyle', ':', 'linewidth',1.2);
text([2550],[4.8], 'Reversal', 'Color', 'k','FontSize', 12, 'fontname','arial')
text([15],[-1.4], 'n=2', 'Color', 'k','FontSize', 12, 'fontname','arial')
ylim([-2 5])
%% CALCULATE CONFIDENCE INTERVAL 95%
xci = [0:0.1:1]; %  Independent Variable
yci = pFA_plateau(WTtouseyoung,:); %  Dependent Variable 
nci = length(yci); % n
ymeanci = nanmean(yci); % Mean Of All Experiments At Each Value Of ?x?
ysdci = nanstd(yci); % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x?
ci95 = tinv([0.025 0.975], nci-1); % Calculate 95% Probability Intervals Of t-Distribution
yci95 = bsxfun(@times, ysdci, ci95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ?x?


mucalc=ymeanci+sqrt(ysdci^2/nci);


figure
plot(xci,ymeanci)                                      % Plot Mean Of All Experiments
hold on
plot(xci, yci95+ymeanci)                                % Plot 95% Confidence Intervals Of All Experiments
hold off
grid


MUTnew=MUTtouseyoung([2, 4:end-2]);
figure
xdist=histogram(pFA_plateau(WTtouseyoung,:),20);

xvals=[0.06:0.06:0.6];

figure
hold all
pd=fitdist(xdist.Values','Binomial');
ydist = pdf(pd,xdist.Data);
plot(xvals,ydist,'LineWidth',2)


histogram(pFA_plateau(MUTnew,:),20);

histogram(pFA_plateau([1 3 20 21],:),20);


figure
hold all
histogram(pFAsli100R(MUTnew,1:maxdtop))

cdf(pFAsli100R(WTtouseyoung,1:maxdtop))


[pv,hv]=ttest2(pFA_plateau(WTtouseyoung,:),pFA_plateau([1 3 20 21],:));



itervalfit=[0.173534:0.338066];
%% FIGURE 3: Fits for individual animals and finds plateau performance
micetoeval=[WTtouseyoung MUTtouseyoung];
maxblock=40;
colors=[[187, 186, 183]; [255, 195, 185]]./255;
plotyes=1; %1 if you want to get figure with individual fits, 0 if not
exampleanimalsR=0; %to plot example Reinforced animals
clear xf yf xfplat yfplat yfplathalf xfplathalf mousename
i=0;
exi=0;
for ii=micetoeval
    i=i+1;
    
    %Reinforced
    xf=[1 xtrialsR(ii,1:maxblock)]'; %x data for fit
    yf=[0 dp_sli100R(ii,1:maxblock)]'; %y data for fit

    %Probe
    xpf=[1 xtrialsP(ii,1:maxblock)]'; %x data for fit
    ypf=[0 dp_slitouseP(ii,1:maxblock)]'; %y data for fit
    
    if sum(~isnan(xf))>=4
        [fitresultR{ii}, gofR{ii}] = mysigmfit(xf, yf);
        %find plateau in reinforced
        clear yplat2 ord closestIndex
        yfit_dp_sli100R{ii}=feval(fitresultR{ii}, xf);
        ord=max(yfit_dp_sli100R{ii})-max(yfit_dp_sli100R{ii})*0.05;
        [minValue,closestIndex]=min(abs(yfit_dp_sli100R{ii}-ord)); %sort(abs(yplat-ordmean))
        if xf(closestIndex)>1000
            xfplat(ii,:)=xf(closestIndex);
            yfplat(ii,:)=yfit_dp_sli100R{ii}(closestIndex);
        else
            vct=xf(~isnan(xf)); vct2=yfit_dp_sli100R{ii}(~isnan(yfit_dp_sli100R{ii}));
            xfplat(ii,:)=vct(end);
            yfplat(ii,:)=vct2(end);
        end
        mousename(ii,:)=char(mouse(ii));

        %find half to plateau
        ordhalf=nanmean(yfit_dp_sli100R{ii});
        [minValuehalf,closestIndexhalf]=min(abs(yfit_dp_sli100R{ii}-ordhalf)); %sort(abs(yplat-ordmean))
        yfplathalf(ii,:)=yfit_dp_sli100R{ii}(closestIndexhalf);
        xfplathalf(ii,:)=xf(closestIndexhalf);        
    end
    if sum(~isnan(xpf))>=4
        %probe
        [fitresultP{ii}, gofP{ii}] = mysigmfit(xpf, ypf);
        yfit_dp_sli100P{ii}=feval(fitresultP{ii}, xpf); 
        
        ordp=max(yfit_dp_sli100P{ii})-max(yfit_dp_sli100P{ii})*0.05;
        [minValuep,closestIndexp]=min(abs(yfit_dp_sli100P{ii}-ordp)); %sort(abs(yplat-ordmean))
        if xpf(closestIndexp)>1000
            xfplatP(ii,:)=xpf(closestIndexp);
            yfplatP(ii,:)=yfit_dp_sli100P{ii}(closestIndexp);
        else
            vctp=xpf(~isnan(xpf)); vct2p=yfit_dp_sli100P{ii}(~isnan(yfit_dp_sli100P{ii}));
            xfplatP(ii,:)=vctp(end);
            yfplatP(ii,:)=vct2p(end);
        end
    else
        yfit_dp_sli100P{ii}=NaN;
    end
    
    if exampleanimalsR==1 && strcmp(mouse(ii),'KF023') || strcmp(mouse(ii),'FZ010') || strcmp(mouse(ii),'AW031') || strcmp(mouse(ii),'AW023')
        gi=figure(10984);
        exi=exi+1;
        subplot(4,4,exi)
        hold all
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
        set(gi, 'position', [1547, 7, 250, 500]);
        h1=plot(fitresultR{ii},'k-',xf(2:end),yf(2:end),'k.');
        set(gca, 'Xtick', [1:4000:8000], 'Xticklabel', [1:4000:8000]);
        xlim([0 10000]); ylim([-1 5]); xlabel('Trials'); ylabel('d'''); title(['(' char(mouse(ii)) ')']); legend off
%         text([300 300],[4.5 4.5], 'Reinforced', 'Color', 'k','FontSize', 12, 'fontname','arial')
%         text([300 300],[3.8 3.8], 'Probe', 'Color', colors(1,:),'FontSize', 12, 'fontname','arial')
        mmi=char(mouse(ii));
        dpexamps.(mmi)=[xf yf yfit_dp_sli100R{ii}];
    end

    if plotyes==1
        gi=figure(10984);
        subplot(4,4,i)
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
        set(gi, 'position', [1547, 7, 553, 736]);
        hold all
        h1=plot(fitresultR{ii},'k-',xf(2:end),yf(2:end),'k.');
            plot(xf(2:end),yf(2:end),'k.');
        h2=plot(fitresultP{ii},'-'); plot(xpf(2:end),ypf(2:end),'.', 'MarkerEdgeColor', colors(1,:));
%         plot(xfplatP(ii,:), yfplatP(ii,:), 'ro') %for plateau
%         plot(xfplathalf(ii,:), yfplathalf(ii,:), 'go') %for halfplateau       
        set(h1, 'linewidth', 2);
            set(h2, {'color'}, num2cell(colors(1,:), 2), 'linewidth', 2);
        set(h2, {'color'}, num2cell(colors(1,:), 2), 'linewidth', 2);
       
        if ismember(ii,WTtouseyoung)
            title(['\color{Blue}' char(mouse(ii))])
        else
            title(char(mouse(ii)))
        end
        
        set(gca, 'Xtick', [0:4000:maxblock*100], 'Xticklabel', [0:4000:maxblock*100]);
        xlim([1 maxblock*100])
        ylim([-1 5])
        xlabel('Trials')
        ylabel('d''')
        legend off
        if i==1
            text([300 300],[4.5 4.5], 'Reinforced', 'Color', 'k','FontSize', 12, 'fontname','arial')
            text([300 300],[3.8 3.8], 'Probe', 'Color', colors(1,:),'FontSize', 12, 'fontname','arial')
        end     
        
%         figure(1012)
%         subplot(1,2,1)
%         hold all
%         h1=plot(fitresultR{ii},'-',xf(2:end),yf(2:end),'.');
%         subplot(1,2,2)
%         hold all
%         h2=plot(fitresultP{ii},'-'); plot(xpf(2:end),ypf(2:end),'.');
    end
end

%% overlap d' reinforced and probe
figure
subplot(1,2,1)
for ii=WTtouseyoung
    xf=[1 xtrialsR(ii,1:maxblock)]'; %x data for fit
    yf=[0 dp_sli100R(ii,1:maxblock)]'; %y data for fit
     h1=plot(fitresultR{ii},'b-',xf(2:end),yf(2:end),'w.');
    hold all
end
for ii=MUTtouseyoung
    xf=[1 xtrialsR(ii,1:maxblock)]'; %x data for fit
    yf=[0 dp_sli100R(ii,1:maxblock)]'; %y data for fit
    h1=plot(fitresultR{ii},'k-',xf(2:end),yf(2:end),'w.');
%     h1=plot(fitresultR{ii},'k-');
    hold all
    title('Reinforced')
    ylabel('d''')
    xlabel('Trials')
end 

subplot(1,2,2)
for ii=WTtouseyoung
    xpf=[1 xtrialsP(ii,1:maxblock)]'; %x data for fit
    ypf=[0 dp_slitouseP(ii,1:maxblock)]'; %y data for fit
     h1=plot(fitresultP{ii},'b-',xpf(2:end),ypf(2:end),'w.');
    hold all
end
for ii=MUTtouseyoung
    xpf=[1 xtrialsP(ii,1:maxblock)]'; %x data for fit
    ypf=[0 dp_slitouseP(ii,1:maxblock)]'; %y data for fit
    h1=plot(fitresultP{ii},'k-',xpf(2:end),ypf(2:end),'w.');
%     h1=plot(fitresultR{ii},'k-');
    hold all
    title('Probe')
    ylabel('d''')
    xlabel('Trials')
    legend off
end 

%% AVERAGE NUMBER OF TRIALS PER DAY
for ii=1:length(mouse)
    trialsperday(ii,:)=round(nanmean(cell2mat(data(ii).ntrials)));
end 


figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, trialsperday(WTtouseyoung,:)', colors(1,:),6,1);
paper_dotbarplot_modif(1, trialsperday(MUTtouseyoung,:)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'WR', 'CA'})
ylabel('Trials/day')
xlim([0 1.5])
pval=mynormandttest(trialsperday(WTtouseyoung), trialsperday(MUTtouseyoung));
title(['p= ' num2str(pval)])


%% GRANT FIGURE PLATEAU FA AND HIT RATE for reinforced
wtmice=[WTtouseyoung];% WTtouseold];
mutmice=[MUTtouseyoung];% MUTtouseold];

clear indxplateau pH_plateau pFA_plateau
% endday=3;
maxnum=200;
for ii=[wtmice mutmice]
    plateau_start=xfplat(ii);
    indxplateau=find(xtrialsR(ii, :)>=plateau_start);
    y_pH=sort(pHsli100R(ii,indxplateau:end));
    y_pFA=sort(pFAsli100R(ii,indxplateau:end));
    y_pH(isnan(y_pH))=[];
    y_pFA(isnan(y_pFA))=[];
    if length(y_pH)>=maxnum
        pH_plateau(ii,:)=nanmean(y_pH(end-(maxnum-1):end));
        pFA_plateau(ii,:)=nanmean(y_pFA(end-(maxnum-1):end));
    else
        pH_plateau(ii,:)=nanmean(y_pH);
        pFA_plateau(ii,:)=nanmean(y_pFA);
    end
end

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, pH_plateau(wtmice,:)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(1, pH_plateau(mutmice,:)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(1.8, pFA_plateau(wtmice,:)', colorsmat2(2,:),6,1);
paper_dotbarplot_modif(2.3, pFA_plateau(mutmice,:)', colorsmat2(2,:),6,1);
set(gca, 'Xtick', [0.5 1 1.8 2.3], 'Xticklabel', {'WR', 'CA'})
ylabel('Action rate')
ylim([0 1.05])
xlim([0 2.8])

% FIGURE STATISTICS
clear tou anov_x anov_y anov_pg anov_g anov_d
tou=[pH_plateau pFA_plateau]; %dp_slitouseP dp_sli100R
anovaWT=wtmice;
anovaMUT=mutmice;
anov_x=[tou(anovaWT,:); tou(anovaMUT,:)];
anov_y=reshape(anov_x,1,[]);
anov_pg=[ones(length(anovaWT),1,1); ones(length(anovaMUT),1,1)*2]';
anov_g=repmat(anov_pg,1,2);
anov_d=[];
for di=3:4
    anov_d=[anov_d ones(1,size(anov_x,1),1)*di];
end

[p,tbl,stats]=anovan(anov_y, {anov_g, anov_d},'varnames',{'group','tone'}, 'model','interaction');
figure; [c, m, h, nms]=multcompare(stats,'Dimension',[1 2]); %, 'ctype', 'bonferroni'

%% GRANT FIGURE PLATEAU FA AND HIT RATE FOR PROBE
clear indxplateau pH_plateauP pFA_plateauP
% endday=3;
maxnum=200;
for ii=[WTtouseyoung MUTtouseyoung]
    plateau_start=xfplat(ii);
    indxplateau=find(xtrialsP(ii, :)>plateau_start);
    y_pH=sort(pHsli100P(ii,indxplateau:end));
    y_pFA=sort(pFAsli100P(ii,indxplateau:end));
    y_pH(isnan(y_pH))=[];
    y_pFA(isnan(y_pFA))=[];
    if length(y_pH)>=maxnum
        pH_plateauP(ii,:)=nanmean(y_pH(end-(maxnum-1):end));
        pFA_plateauP(ii,:)=nanmean(y_pFA(end-(maxnum-1):end));
    else
        pH_plateauP(ii,:)=nanmean(y_pH);
        pFA_plateauP(ii,:)=nanmean(y_pFA);
    end
end

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, pH_plateauP(WTtouseyoung,:)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(1, pH_plateauP(MUTtouseyoung,:)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(1.8, pFA_plateauP(WTtouseyoung,:)', colorsmat2(2,:),6,1);
paper_dotbarplot_modif(2.3, pFA_plateauP(MUTtouseyoung,:)', colorsmat2(2,:),6,1);
set(gca, 'Xtick', [0.5 1 1.8 2.3], 'Xticklabel', {'WR', 'CA'})
ylabel('Action rate')
ylim([0 1.05])
xlim([0 2.8])

% FIGURE STATISTICS
clear tou anov_x anov_y anov_pg anov_g anov_d
tou=[pH_plateauP pFA_plateauP]; %dp_slitouseP dp_sli100R
anovaWT=WTtouseyoung;
anovaMUT=MUTtouseyoung;
anov_x=[tou(anovaWT,:); tou(anovaMUT,:)];
anov_y=reshape(anov_x,1,[]);
anov_pg=[ones(length(anovaWT),1,1); ones(length(anovaMUT),1,1)*2]';
anov_g=repmat(anov_pg,1,2);
anov_d=[];
for di=3:4
    anov_d=[anov_d ones(1,size(anov_x,1),1)*di];
end

[p,tbl,stats]=anovan(anov_y, {anov_g, anov_d},'varnames',{'group','tone'}, 'model','interaction');
figure; [c, m, h, nms]=multcompare(stats,'Dimension',[1 2]); %, 'ctype', 'bonferroni'

%% FIGURE GRANT 1: individual pH pFA action rates for reinforced and probe
clear xp yp semp
maxdtop=40; 
colorsmat2=[96, 211, 148; 242, 143, 59]./255; % 238, 96, 85 . [[155, 204, 85]./255; [228, 87, 46]./255; [179, 201, 174]./255; [239, 149, 137]./255];
exampleanimalsR=0;
gi=figure(1029);
set(gi, 'position', [1547, 7, 581, 736]); 
% set(gcf, 'Position',  [100, 100, 2000, 300])
mitoplot=[WTtouseyoung MUTtouseyoung]; %[WTtouseold MUTtouseold]; %new2MUTtouse; %[WTtouseyoung WTtouseold MUTtouseyoung MUTtouseold];
si=0;
exi=0;
for ii=mitoplot %[mitoplot] %MUTtouseyoung AllMUTold WTtouseyoung
    si=si+1; %Reinforced
    
    figure(1028)
    subplot(4,4,si); hold all
    set(gca,'TickDir','out','FontSize',12,'fontname','arial');
    if ismember(ii,WTtouseyoung)
        title(['\color{Blue}' char(mouse(ii))])
    else
        title(char(mouse(ii)))
    end
    ylabel('Action rate'); xlabel('Reinforced trials');
    plot(xtrialsR(ii,1:maxdtop), smooth(pHsli100R(ii,1:maxdtop)),'-','color', colorsmat2(1,:), 'linewidth', 2)
    plot(xtrialsR(ii,1:maxdtop), smooth(pFAsli100R(ii,1:maxdtop)),'-','color', colorsmat2(2,:), 'linewidth', 2)
%     text([200 200],[0.2 0.2], 'Target', 'Color', colorsmat2(1,:),'FontSize', 12, 'fontname','arial')
%     text([200 200], [0.08 0.08], 'Foil', 'Color', colorsmat2(2,:),'FontSize', 12, 'fontname','arial')
    set(gca, 'Xtick', [1 (maxdtop*100)/2 maxdtop*100], 'Xticklabel', [1 (maxdtop*100)/2  maxdtop*100]);
    xlim([1 maxdtop*100]); ylim([0 1.02]);
    line([0 10000], [0.98 0.98], 'color', 'k');
     if exampleanimalsR==1 && strcmp(mouse(ii),'KF023') || strcmp(mouse(ii),'AM024') || strcmp(mouse(ii),'AW031') || strcmp(mouse(ii),'AW023')
        gi=figure(109898);
        exi=exi+1;
        subplot(2,2,exi)
        hold all
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
%         set(gi, 'position', [1547, 7, 250, 500]);
        ylabel('Action rate'); xlabel('Reinforced trials'); 
        plot(xtrialsR(ii,1:maxdtop),smooth(pHsli100R(ii,1:maxdtop)),'-','color', colorsmat2(1,:), 'linewidth', 2)
        plot(xtrialsR(ii,1:maxdtop),smooth(pFAsli100R(ii,1:maxdtop)),'-','color', colorsmat2(2,:), 'linewidth', 2)
        text([200],[0.2], 'Target', 'Color', colorsmat2(1,:),'FontSize', 12, 'fontname','arial')
        text([200], [0.08], 'Foil', 'Color', colorsmat2(2,:),'FontSize', 12, 'fontname','arial')
        set(gca, 'Xtick', [0:4000:12000], 'Xticklabel', [0:4000:12000]);
        xlim([1 maxdtop*100]); ylim([0 1.02]); title(['(' char(mouse(ii)) ')']); legend off
        mmi=char(mouse(ii));
        arexamps.(mmi)=[xtrialsR(ii,1:maxdtop)' pHsli100R(ii,1:maxdtop)' pFAsli100R(ii,1:maxdtop)'];
    end
% %     
    %Probe
%     figure(1031)
%     subplot(4,4,si); hold all
%     set(gca,'TickDir','out','FontSize',12,'fontname','arial');
%     if ismember(ii,WTtouseyoung)
%         title(['\color{Blue}' char(mouse(ii))])
%     else
%         title(char(mouse(ii)))
%     end
%     ylabel('Action rate'); xlabel('Probe trials'); xlim([1 8000]); ylim([0 1.02]);
%     plot(xtrialsP(ii,1:maxdtop),smooth(pHsli100P(ii,1:maxdtop)),'-','color', colorsmat2(1,:), 'linewidth', 2)
%     plot(xtrialsP(ii,1:maxdtop),smooth(pFAsli100P(ii,1:maxdtop)),'-','color', colorsmat2(2,:), 'linewidth', 2)
% %     text([200 200],[0.2 0.2], 'Target', 'Color', colorsmat2(1,:),'FontSize', 12, 'fontname','arial')
% %     text([200 200], [0.08 0.08], 'Foil', 'Color', colorsmat2(2,:),'FontSize', 12, 'fontname','arial')
%     set(gca, 'Xtick', [0:4000:8000], 'Xticklabel', [0:4000:8000]);
end

%% overlay action rates 
figure
for ii=WTtouseyoung
    subplot(2,2,1)
    hold all
    plot(xtrialsR(ii,1:maxdtop), smooth(pHsli100R(ii,1:maxdtop)),'-','color', colorsmat2(1,:), 'linewidth', 2)
    plot(xtrialsR(ii,1:maxdtop), smooth(pFAsli100R(ii,1:maxdtop)),'-','color', colorsmat2(2,:), 'linewidth', 2)
    title('Reinforced'); ylim([0 1]);
    
    subplot(2,2,2)
    hold all
    plot(xtrialsP(ii,1:maxdtop),smooth(pHsli100P(ii,1:maxdtop)),'-','color', colorsmat2(1,:), 'linewidth', 2)
    plot(xtrialsP(ii,1:maxdtop),smooth(pFAsli100P(ii,1:maxdtop)),'-','color', colorsmat2(2,:), 'linewidth', 2)
    title('Probe'); ylim([0 1]);
end

for ii=MUTtouseyoung
    subplot(2,2,3)
    hold all
    plot(xtrialsR(ii,1:maxdtop), smooth(pHsli100R(ii,1:maxdtop)),'-','color', colorsmat2(1,:), 'linewidth', 2)
    plot(xtrialsR(ii,1:maxdtop), smooth(pFAsli100R(ii,1:maxdtop)),'-','color', colorsmat2(2,:), 'linewidth', 2)
    title('Reinforced'); ylim([0 1]); ylabel('Action Rates'); xlabel('Trial');
    
    subplot(2,2,4)
    hold all
    plot(xtrialsP(ii,1:maxdtop),smooth(pHsli100P(ii,1:maxdtop)),'-','color', colorsmat2(1,:), 'linewidth', 2)
    plot(xtrialsP(ii,1:maxdtop),smooth(pFAsli100P(ii,1:maxdtop)),'-','color', colorsmat2(2,:), 'linewidth', 2)
    title('Probe'); ylim([0 1]);
end


%% Learning hit rate for old animals
clear HitRthr FARthr
hitthresh=0.92; %for young =0.98
fathresh=0.82;  %for young =0.92

for ii=[WTtouseyoung MUTtouseyoung]
    clear findthrhit hitvar_pHsli100R hitvar favar_pHsli100R favar findthrfa
    hitvar_pHsli100R=pHsli100R(ii,:);
    favar_pHsli100R=pFAsli100R(ii,:);
    hitvar=hitvar_pHsli100R(~isnan(hitvar_pHsli100R));
    favar=favar_pHsli100R(~isnan(favar_pHsli100R));
    findthrhit=find(smooth(hitvar)>=hitthresh);
    findthrfa=find(smooth(favar)>=fathresh);
    if ~isempty(findthrhit)
        HitRthr(ii,:)=xtrialsR(ii,findthrhit(1));
    else
        HitRthr(ii,:)=NaN;
    end
    
    if ~isempty(findthrfa)
        FARthr(ii,:)=xtrialsR(ii,findthrfa(1));
    else
        FARthr(ii,:)=NaN;
    end
end

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, HitRthr(WTtouseyoung,:)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(1, HitRthr(MUTtouseyoung,:)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(1.8, FARthr(WTtouseyoung,:)', colorsmat2(2,:),6,1);
paper_dotbarplot_modif(2.3, FARthr(MUTtouseyoung,:)', colorsmat2(2,:),6,1);
set(gca, 'Xtick', [0.5 1 1.8 2.3], 'Xticklabel', {'WR', 'CA'})
ylabel('Max threshold trial #')
xlim([0 2.8])
% [pval1]=mynormandttest(HitRthr(WTtouseold,:), HitRthr(MUTtouseold,:));
% title(['p=' num2str(pval1)])
% ylim([0 2.8])

% STATISTICS
stats_x=[HitRthr(WTtouseyoung,:); HitRthr(MUTtouseyoung,:); FARthr(WTtouseyoung,:); FARthr(MUTtouseyoung,:)]';
stats_g=[ones(1,length(WTtouseyoung)), ones(1,length(MUTtouseyoung))*2, ones(1,length(WTtouseyoung))*3, ones(1,length(MUTtouseyoung))*4];
stats_cond=[repmat({'H'}, 1, length(WTtouseyoung)), repmat({'H'}, 1, length(MUTtouseyoung)), repmat({'FA'}, 1, length(WTtouseyoung)), repmat({'FA'}, 1, length(MUTtouseyoung))];
kstest(stats_x);

[p,tbl,stats]=anovan(stats_x, {stats_g, stats_cond},'varnames',{'group','action'}, 'model','interaction');
figure
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'ctype', 'bonferroni')


[p,tbl,stats] = kruskalwallis(stats_x, stats_g)

%% NOT IN USE: To bin probe in additional bins

dp_slibintouseP=NaN(32,120);
dp_WTslibintouseP=NaN(20,120);
dp_MUTslibintouseP=NaN(20,120);
winincrease=600;
ini=0;
for sc=1:120
    xbin(:,sc)=ini+winincrease;
    fivectwt=find(xtrialsP(WTtouseyoung,:)>ini & xtrialsP(WTtouseyoung,:)<=ini+winincrease);
    fivectmut=find(xtrialsP(MUTtouseyoung,:)>ini & xtrialsP(MUTtouseyoung,:)<=ini+winincrease);
    if ~isempty(fivectwt)
        dp_WTslibintouseP(1:length(fivectwt),sc)=dp_slitouseP(fivectwt);
    end
    if ~isempty(fivectmut)
        dp_MUTslibintouseP(1:length(fivectmut),sc)=dp_slitouseP(fivectmut);
    end 
    ini=ini+winincrease;
end
% figure
% hold all
% plot(xbin(1:40), dp_WTslibintouseP(:,1:40), 'ko')
% plot(xbin(1:40), dp_MUTslibintouseP(:,1:40), 'ro')
%SEM
% clear xr xpwt xpmut 
% xr=dp_sli100R;
% xpwt=dp_WTslibintouseP;
% xpmut=dp_MUTslibintouseP;
% 
% clear SEMRY_WT SEMRY_MUT SEMPY_WT SEMPY_MUT
% for ai=1:120
%     wtsR=sum(~isnan(xr(WTtouseyoung,ai)));
%     mutsR=sum(~isnan(xr(MUTtouseyoung,ai)));
%     wtsP=sum(~isnan(xpwt(:,ai)));
%     mutsP=sum(~isnan(xpmut(:,ai)));
%     SEMRY_WT(ai,:)=nanstd(xr(WTtouseyoung,ai))/sqrt(wtsR);
%     SEMRY_MUT(ai,:)=nanstd(xr(MUTtouseyoung,ai))/sqrt(mutsR);
%     SEMPY_WT(ai,:)=nanstd(xpwt(:,ai))/sqrt(wtsP);
%     SEMPY_MUT(ai,:)=nanstd(xpmut(:,ai))/sqrt(mutsP);
% end
%% FIGURE 1: d primes per trial GOOD
wtmice=[WTtouseyoung]; % WTtouseyoung]; %WTtouseyoung; %WTtouseold;
mutmice=[MUTtouseyoung]; % MUTtouseold]; %MUTtouseyoung; %MUTtouseold;

%SEM
clear xr xpwt xpmut 
xr=dp_sli100R;
xp=dp_slitouseP;

clear SEMRY_WT SEMRY_MUT SEMPY_WT SEMPY_MUT
for ai=1:120
    wtsR=sum(~isnan(xr(wtmice,ai)));
    mutsR=sum(~isnan(xr(mutmice,ai)));
    wtsP=sum(~isnan(xp(wtmice,ai)));
    mutsP=sum(~isnan(xp(mutmice,ai)));
    SEMRY_WT(ai,:)=nanstd(xr(wtmice,ai))/sqrt(wtsR);
    SEMRY_MUT(ai,:)=nanstd(xr(mutmice,ai))/sqrt(mutsR);
    SEMPY_WT(ai,:)=nanstd(xp(:,ai))/sqrt(wtsP);
    SEMPY_MUT(ai,:)=nanstd(xp(:,ai))/sqrt(mutsP);
end

maxdtop=35;
maxdnext=30; %38;
extraday=30;
jumps=1;
%reinf
xtoplotRwt=[nanmean(xtrialsR(wtmice, 1:maxdtop))];
ytoplotRwt=[nanmean(dp_sli100R(wtmice, 1:maxdtop))];
xtoplotRmut=[nanmean(xtrialsR(mutmice, 1:maxdtop))];
ytoplotRmut=[nanmean(dp_sli100R(mutmice, 1:maxdtop))];
%probe
xtoplotPwt=[nanmean(xtrialsP(wtmice, 1:maxdtop))];
ytoplotPwt=[nanmean(dp_slitouseP(wtmice, 1:maxdtop))];
xtoplotPmut=[nanmean(xtrialsP(mutmice, 1:maxdtop))];
ytoplotPmut=[nanmean(dp_slitouseP(mutmice,1:maxdtop))];

figure
% subplot(1,3,1)
% set(gca,'TickDir','out','FontSize',12,'fontname','arial');
% hold all
% plot(xtoplotRwt, ytoplotRwt, 'b.'); 
% plot(xtoplotRmut, ytoplotRmut, 'k.'); 
% plot(xtoplotPwt, ytoplotPwt, '.', 'color', colors(1, :)); 
% plot(xtoplotPmut, ytoplotPmut, '.', 'color', colors(2, :)); 
[fitresultWT, gofWT] = mysigmfit([1 xtoplotRwt]', [0 ytoplotRwt]');
[fitresultMUT, gofMUT] = mysigmfit([1 xtoplotRmut]', [0 ytoplotRmut]');
[fitresultWTP, gofWTP] = mysigmfit([1 xtoplotPwt]',[0 ytoplotPwt]');
[fitresultMUTP, gofMUP] = mysigmfit([1 xtoplotPmut]', [0 ytoplotPmut]');
% a=plot(fitresultWT,'b-'); set(a,'linewidth', 2);
% a1=plot(fitresultMUT,'k-'); set(a1,'linewidth', 2);
% b=plot(fitresultWTP,'b-'); set(b,{'color'},num2cell(colors(1, :),2),'linewidth', 2);
% b1=plot(fitresultMUTP,'k-'); set(b1,{'color'},num2cell(colors(2, :),2),'linewidth', 2);
% text([300 300],[3.7 3.7], 'WR', 'Color', 'b','FontSize', 12, 'fontname','arial')
% text([300 300],[3.3 3.3], 'CA', 'Color', 'k','FontSize', 12, 'fontname','arial')
% text([300 300],[2.8 2.8], 'WR', 'Color', 'b','FontSize', 12, 'fontname','arial')
% text([300 300],[2.4 2.4], 'CA', 'Color', 'k','FontSize', 12, 'fontname','arial')
% ylim([-1 5])
% xlim([1 maxdtop*100])
% ylabel('d''')
% xlabel('Trials')
% legend off
subplot(1,2,1)
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
line([5000 5000],[-1 5],'color', 'k', 'linestyle', ':')
errorbar(xtoplotRwt(1:jumps:maxdnext), ytoplotRwt(1:jumps:maxdnext), SEMRY_WT(1:jumps:maxdnext)', 'k.', 'LineStyle','none');
errorbar(xtoplotRmut(1:jumps:maxdnext), ytoplotRmut(1:jumps:maxdnext), SEMRY_MUT(1:jumps:maxdnext)', '.', 'color', plotscolor, 'LineStyle','none');
errorbar(xtoplotRwt(extraday), ytoplotRwt(extraday), SEMRY_WT(extraday)', 'k.', 'LineStyle','none');
errorbar(xtoplotRmut(extraday), ytoplotRmut(extraday), SEMRY_MUT(extraday)', '.', 'color', plotscolor, 'LineStyle','none');
set(gca, 'Xtick', [0:2000:8000], 'Xticklabel', [0:2000:8000]);
a=plot(fitresultWT,'k-'); set(a,'linewidth', 2);
a1=plot(fitresultMUT,'-'); set(a1,'color',plotscolor,'linewidth', 2);
text([300 300],[4.8 4.8], 'saline', 'Color', 'k','FontSize', 12, 'fontname','arial')
text([300 300],[4.5 4.5], 'dsp4', 'Color', plotscolor,'FontSize', 12, 'fontname','arial')
ylim([-1 5])
legend off
xlim([1 maxdtop*100])
ylabel('d''')
xlabel('Trials')
title('Reinforced')
subplot(1,2,2)
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
line([5000 5000],[-1 4],'color', 'k', 'linestyle', ':')
errorbar(xtoplotPwt(1:9), ytoplotPwt(1:9), SEMPY_WT(1:9)','.', 'color', 'k', 'LineStyle','none');
errorbar(xtoplotPmut(1:7), ytoplotPmut(1:7), SEMPY_MUT(1:7)','.', 'color',plotscolor, 'LineStyle','none');
errorbar(xtoplotPwt(13), ytoplotPwt(13), SEMPY_WT(13)','.', 'color','k', 'LineStyle','none');
errorbar(xtoplotPmut(11), ytoplotPmut(11), SEMPY_MUT(11)','.', 'color',plotscolor, 'LineStyle','none');
b=plot(fitresultWTP,'k-'); set(b,'color','k','linewidth', 2);
b1=plot(fitresultMUTP,'-'); set(b1,'color',plotscolor,'linewidth', 2);
text([300 300],[4.8 4.8], 'saline', 'Color', 'k','FontSize', 12, 'fontname','arial')
text([300 300],[4.5 4.5], 'dsp4', 'Color', plotscolor, 'FontSize', 12, 'fontname','arial')
ylim([-1 5])
xlim([1 maxdtop*100])
set(gca, 'Xtick', [0:2000:8000], 'Xticklabel', [0:2000:8000]);
ylabel('d''')
xlabel('Trials')
title('Probe')
legend off



%% FIGURE 1: d primes per trial ONLY REINFORCED For grant
colors=[[187, 186, 183]; [255, 195, 185]]./255;

clear xtoplotRwt ytoplotRwt xtoplotRmut ytoplotRmut fitresultWT fitresultMUT fitresultWTo fitresultMUTo
wtmice=WTtouseyoung; %WTtouseold; %WTtouseyoung; %WTtouseold;
mutmice=MUTtouseyoung; %MUTtouseold; %MUTtouseyoung; %MUTtouseold;
% wtmiceo=WTtouseold;
% mutmiceo=MUTtouseold;

%SEM
clear xr xpwt xpmut 
xr=dp_slitouseP; %dp_sli100R;
xx=xtrialsP;% xtrialsR;

clear SEMRY_WT SEMRY_MUT SEMPY_WT SEMPY_MUT
for ai=1:80
    wtsR=sum(~isnan(xr(wtmice,ai)));
    mutsR=sum(~isnan(xr(mutmice,ai)));
    SEMRY_WT(ai,:)=nanstd(xr(wtmice,ai))/sqrt(wtsR);
    SEMRY_MUT(ai,:)=nanstd(xr(mutmice,ai))/sqrt(mutsR);
    
%     wtsRo=sum(~isnan(xr(wtmiceo,ai)));
%     mutsRo=sum(~isnan(xr(mutmiceo,ai)));
%     SEMRY_WTo(ai,:)=nanstd(xr(wtmiceo,ai))/sqrt(wtsRo);
%     SEMRY_MUTo(ai,:)=nanstd(xr(mutmiceo,ai))/sqrt(mutsRo);
end

%for younger cohorot
maxdtop=80;
maxdnext=38;
extraday=53;
jumps=1;
xtoplotRwt=[nanmean(xx(wtmice, 1:maxdtop))];
ytoplotRwt=[nanmean(xr(wtmice, 1:maxdtop))];
xtoplotRmut=[nanmean(xx(mutmice, 1:maxdtop))];
ytoplotRmut=[nanmean(xr(mutmice, 1:maxdtop))];
[fitresultWT, gofWT] = mysigmfit([1 xtoplotRwt]', [0 ytoplotRwt]');
[fitresultMUT, gofMUT] = mysigmfit([1 xtoplotRmut]', [0 ytoplotRmut]');

%older cohort
% maxdtopo=120;
% maxdnexto=38;
% extradayo=73;
% xtoplotRwto=[nanmean(xx(wtmiceo, 1:maxdtopo))];
% ytoplotRwto=[nanmean(xr(wtmiceo, 1:maxdtopo))];
% xtoplotRmuto=[nanmean(xx(mutmiceo, 1:maxdtopo))];
% ytoplotRmuto=[nanmean(xr(mutmiceo, 1:maxdtopo))];
% [fitresultWTo, gofWTo] = mysigmfit([1 xtoplotRwto]', [0 ytoplotRwto]');
% [fitresultMUTo, gofMUTo] = mysigmfit([1 xtoplotRmuto]', [0 ytoplotRmuto]');


%%
figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
line([5000 5000],[-1 4],'color', 'k', 'linestyle', ':')
errorbar(xtoplotRwt(1:jumps:maxdnext), ytoplotRwt(1:jumps:maxdnext), SEMRY_WT(1:jumps:maxdnext)', 'k.', 'LineStyle','none');
errorbar(xtoplotRmut(1:jumps:maxdnext), ytoplotRmut(1:jumps:maxdnext), SEMRY_MUT(1:jumps:maxdnext)', 'r.', 'LineStyle','none');
errorbar(xtoplotRwt(extraday), ytoplotRwt(extraday), SEMRY_WT(extraday)', 'k.', 'LineStyle','none');
errorbar(xtoplotRmut(extraday), ytoplotRmut(extraday), SEMRY_MUT(extraday)', 'r.', 'LineStyle','none');
a=plot(fitresultWT,'k-'); set(a,'linewidth', 2);
a1=plot(fitresultMUT,'r-'); set(a1,'linewidth', 2);

errorbar(xtoplotRwto(1:jumps:maxdnexto), ytoplotRwto(1:jumps:maxdnexto), SEMRY_WTo(1:jumps:maxdnexto)', 'ko', 'LineStyle','none', 'Markersize', 3, 'Markerfacecolor', 'w');
errorbar(xtoplotRmuto(1:jumps:maxdnexto), ytoplotRmuto(1:jumps:maxdnexto), SEMRY_MUTo(1:jumps:maxdnexto)', 'ro', 'LineStyle','none', 'Markersize', 3, 'Markerfacecolor', 'w');
errorbar(xtoplotRwto(extraday), ytoplotRwto(extraday), SEMRY_WTo(extraday)', 'ko', 'LineStyle','none','Markersize', 3, 'Markerfacecolor', 'w');
errorbar(xtoplotRmuto(extraday), ytoplotRmuto(extraday), SEMRY_MUTo(extraday)', 'ro', 'LineStyle','none','Markersize', 3, 'Markerfacecolor', 'w');
errorbar(xtoplotRwto(extradayo), ytoplotRwto(extradayo), SEMRY_WTo(extradayo)', 'ko','LineStyle','none','Markersize', 3, 'Markerfacecolor', 'w');
errorbar(xtoplotRmuto(extradayo), ytoplotRmuto(extradayo), SEMRY_MUTo(extradayo)', 'ro','LineStyle','none','Markersize', 3, 'Markerfacecolor', 'w');
ao=plot(fitresultWTo,'k--'); set(ao,'linewidth', 2);
a1o=plot(fitresultMUTo,'r--');set(a1o,'linewidth', 2);
% legend({'', '',  'Wt (6-8mo)', 'APP^+ (6-8mo)','Wt (10-12mo)','APP^+ (10-12mo)'})
% text([300 300],[3.9 3.9], 'Wt (6-8mo)', 'Color', 'k','FontSize', 12, 'fontname','arial')
% text([300 300],[3.5 3.5], 'APP^+ (6-8mo)', 'Color', 'r','FontSize', 12, 'fontname','arial')
% text([300 300],[3.1 3.1], 'Wt (10-12mo)', 'Color', colors(1, :),'FontSize', 12, 'fontname','arial')
% text([300 300],[2.7 2.7], 'APP^+ (10-12mo)', 'Color', colors(2, :),'FontSize', 12, 'fontname','arial')
ylim([-1 4]); xlim([1 6000]); ylabel('d'''); xlabel('Trials'); title('Probe'); %legend off
%% Quantification d' using only first 10 trials of the total 20 probe
normRdpWT=NaN(10,30);
normRdpMUT=NaN(11,30);
xfWT=NaN(10,30);
xfMUT=NaN(11,30);

normRdpWTall=NaN(10,30);
normRdpMUTall=NaN(11,30);
xfWTall=NaN(10,30);
xfMUTall=NaN(11,30);

figure
subplot(1,2,1)
hold all
i=0;
for ii=[WTtouseyoung]
    i=i+1;
    plot(ppi10_xtrialP{ii},(ppi10_dp_P{ii}*100)/2.5631, 'ko-');
    normRdpWT(i,1:length((ppi10_dp_P{ii}*100)/2.5631))=(ppi10_dp_P{ii}*100)/2.5631;
    xfWT(i,1:length((ppi10_dp_P{ii}*100)/2.5631))=ppi10_xtrialP{ii};
end 
i=0;
for ii=[MUTtouseyoung]
    i=i+1;
    plot(ppi10_xtrialP{ii},(ppi10_dp_P{ii}*100)/2.5631, 'ro-');
    normRdpMUT(i,1:length((ppi10_dp_P{ii}*100)/2.5631))=(ppi10_dp_P{ii}*100)/2.5631;
    xfMUT(i,1:length((ppi10_dp_P{ii}*100)/2.5631))=ppi10_xtrialP{ii};
end 
plot(nanmean(xfWT), nanmean(normRdpWT), 'k', 'linewidth', 3)
plot(nanmean(xfMUT), nanmean(normRdpMUT), 'r', 'linewidth', 3)
ylim([-100 105])
ylabel('Normalized d''')
xlabel('Trials')
title('First 10 Probe trials')
xlim([0 6000])

subplot(1,2,2)
hold all
i=0;
for ii=[WTtouseyoung]
    i=i+1;
    plot(xtrialsP(ii,:),(dp_slitouseP(ii,:)*100)/3.2897, 'ko-');
    normRdpWTall(i,1:length(dp_slitouseP(ii,:)))=(dp_slitouseP(ii,:)*100)/3.2897;
    xfWTall(i,1:length(dp_slitouseP(ii,:)))=xtrialsP(ii,:);
end 
i=0;
for ii=[MUTtouseyoung]
    i=i+1;
    plot(xtrialsP(ii,:),(dp_slitouseP(ii,:)*100)/3.2897, 'ro-');
    normRdpMUTall(i,1:length(dp_slitouseP(ii,:)))=(dp_slitouseP(ii,:)*100)/3.2897;
    xfMUTall(i,1:length(dp_slitouseP(ii,:)))=xtrialsP(ii,:);
end 
plot(nanmean(xfWTall), nanmean(normRdpWTall), 'k', 'linewidth', 3)
plot(nanmean(xfMUTall), nanmean(normRdpMUTall), 'r', 'linewidth', 3)
ylim([-100 105])
ylabel('d''')
xlabel('Trials')
title('20 Probe trials')
xlim([0 6000])




%%


figure
hold all

% ylim([-100 105])
ylabel('Normalized d''')
xlabel('Trials')
title('First 10 Probe trials')
xlim([0 6000])


figure
hold all
plot(nanmean(xfWT), nanmean(normRdpWT), 'k')
plot(nanmean(xfMUT), nanmean(normRdpMUT), 'r')
% ylim([-100 105])
ylabel('Normalized d''')
xlabel('Trials')
title('First 10 Probe trials')
xlim([0 6000])

%% FIGURE 1 STATISTICS
clear tou anov_x anov_y anov_pg anov_g anov_d
mindtop=1;
maxdtop=6; %13 for probe
tou=dp_slitouseP(:,mindtop:maxdtop); %dp_slitouseP   dp_sli100R
anovaWT=[WTtouseyoung WTtouseold];
anovaMUT=[MUTtouseyoung MUTtouseold];
anov_x=[tou(WTtouseyoung, :); tou(WTtouseold,:); tou(MUTtouseyoung,:); tou(MUTtouseold,:)];
anov_y=reshape(anov_x,1,[]);
anov_pg=[ones(length(WTtouseyoung),1,1); ones(length(WTtouseold),1,1)*2; ones(length(MUTtouseyoung),1,1)*3; ones(length(MUTtouseold),1,1)*4]';
anov_g=repmat(anov_pg,1,(maxdtop-mindtop)+1);
anov_d=[];
for di=1:(maxdtop-mindtop)+1
    anov_d=[anov_d ones(1,size(anov_x,1),1)*di];
end

anov_g(isnan(anov_y))=[];
anov_d(isnan(anov_y))=[];
anov_y(isnan(anov_y))=[];

[p,tbl,stats]=anovan(anov_y, {anov_g, anov_d},'varnames',{'group','day'}, 'model','interaction');
figure
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'ctype', 'bonferroni')

%%

[p,tbl,stats]=kruskalwallis(anov_y, anov_d);
[c, m, h, nms]=multcompare(stats, 'alpha', 0.05, 'ctype', 'bonferroni')


%% FIGURE 1 STATISTICS
wtstat=dp_sli100R(WTtouseyoung, 1:maxdtop);
mutstat=dp_sli100R(MUTtouseyoung, 1:maxdtop);
datastat=[wtstat; mutstat];
anistats=[repmat({'Wt'}, length(WTtouseyoung),1); repmat({'Mut'},length(MUTtouseyoung),1)];
anistatsnum=[ones(length(WTtouseyoung),1); ones(length(MUTtouseyoung),1)*2];

t=table(anistats,datastat(:,1),datastat(:,2),datastat(:,3),datastat(:,4),datastat(:,5),datastat(:,6),datastat(:,7),datastat(:,8),...
    datastat(:,9),datastat(:,10),datastat(:,11),datastat(:,12),datastat(:,13),datastat(:,14),datastat(:,15),datastat(:,16),datastat(:,17),...
    datastat(:,18),datastat(:,19),datastat(:,20),datastat(:,21),datastat(:,22),datastat(:,23),datastat(:,24),datastat(:,25),datastat(:,26),datastat(:,27),...
    datastat(:,28),datastat(:,29),datastat(:,30),'Variablenames', {'animals','meas1','meas2','meas3','meas4','meas5','meas6', 'meas7', 'meas8', 'meas9', 'meas10',...
    'meas11', 'meas12', 'meas13', 'meas14', 'meas15', 'meas16', 'meas17', 'meas18', 'meas19', 'meas20', 'meas21', 'meas22', 'meas23', 'meas24', 'meas25',...
    'meas26', 'meas27', 'meas28', 'meas29', 'meas30'});

Days=[1:30]';
rm=fitrm(t,'meas1-meas30 ~ animals','WithinDesign',Days,'WithinModel','orthogonalcontrasts');
ranovatbl=ranova(rm);

%second part
t2=table(anistats,datastat(:,31),datastat(:,32),datastat(:,33),datastat(:,34),datastat(:,35),datastat(:,36),datastat(:,37),datastat(:,38),...
    datastat(:,39),datastat(:,40),datastat(:,41),datastat(:,42),datastat(:,43),datastat(:,44),datastat(:,45),datastat(:,46),datastat(:,47),...
    datastat(:,48),datastat(:,49),datastat(:,50),datastat(:,51),datastat(:,52),datastat(:,53),datastat(:,54),datastat(:,55),datastat(:,56),datastat(:,57),...
    datastat(:,58),datastat(:,59),datastat(:,60),'Variablenames', {'animals','meas31','meas32','meas33','meas34','meas35','meas36', 'meas37', 'meas38', 'meas39', 'meas40',...
    'meas41', 'meas42', 'meas43', 'meas44', 'meas45', 'meas46', 'meas47', 'meas48', 'meas49', 'meas50', 'meas51', 'meas52', 'meas53', 'meas54', 'meas55',...
    'meas56', 'meas57', 'meas58', 'meas59', 'meas60'});

Days2=[31:60]';
rm2=fitrm(t2,'meas31-meas60 ~ animals','WithinDesign',Days2); %,'WithinModel','orthogonalcontrasts'
ranovatbl2=ranova(rm2);
mult2=multcompare(rm2, 'animals', 'By', 'meas31-meas60','ComparisonType','bonferroni');




%% SLOPE of dp fit as a measure of learning 
colors2=[0, 184, 245; 254, 221, 11]./255;

clear fitslopewt fitslopemut
i=0;
for ii=WTtouseyoung
    i=i+1;
    fitslopewt(i,:)=fitresultR{ii}.c/fitresultR{ii}.d;
end
i=0;
for ii=MUTtouseyoung
    i=i+1;
    fitslopemut(i,:)=fitresultR{ii}.c/fitresultR{ii}.d;
end

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5,fitslopewt', colors2(1,:),6,1);
paper_dotbarplot_modif(1, fitslopemut', colors2(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'WR', 'CA'})
ylabel('Slope of sigmoid fit')
[pval]=mynormandttest(fitslopewt,fitslopemut);
title([' p= ' num2str(pval)])
% ylim([1000 10000])
xlim([0 1.5])



%% save vector of  and halfplateau DIFFERENCE IN PLATEAU trial #
colors2(1,:)=[0,0,0];
colors2(2,:)=plotscolor;

clear plateauvect
% [pval]=mynormandttest(xfplat(WTtouseyoung),xfplat(MUTtouseyoung));
% [pval]=mynormandttest(xfplat(WTtouseold),xfplat(MUTtouseold));
platvect_good.mouse=mousename;
platvect_good.plateau=xfplat;
platvect_good.halfplateau=xfplathalf;

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, xfplat([WTtouseyoung])', colors2(1,:),6,1);
paper_dotbarplot_modif(1, xfplat([MUTtouseyoung])', colors2(2,:),6,1);
paper_dotbarplot_modif(1.8, xfplatP([WTtouseyoung])', colors2(1,:),6,1);
paper_dotbarplot_modif(2.3, xfplatP([MUTtouseyoung])', colors2(2,:),6,1);
set(gca, 'Xtick', [0.5 1 1.8 2.3], 'Xticklabel', {'WR', 'CA'})
ylabel('Plateau Trial #')
% [pval]=mynormandttest(xfplat([WTtouseyoung WTtouseold]),xfplat([MUTtouseyoung MUTtouseold]));
% title([' p= ' num2str(pval)])
% ylim([1000 11000])
xlim([0 2.8])


%only 6-8mo animals HALF PLATEAU
figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, xfplathalf([WTtouseyoung])', colors2(1,:),6,1);
paper_dotbarplot_modif(1, xfplathalf([MUTtouseyoung])', colors2(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'})
ylabel('Half-plateau Trial #')
% [pval]=mynormandttest(xfplathalf([WTtouseyoung]),xfplathalf([MUTtouseyoung]));
% title([' p= ' num2str(pval)])
% ylim([1000 10000])
xlim([0 1.5])


% substraction from plateau and halfplateau 6-8mo 

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, xfplat([WTtouseyoung])'-xfplathalf([WTtouseyoung])', 'k',6,1);
paper_dotbarplot_modif(1, xfplat([MUTtouseyoung])'-xfplathalf([MUTtouseyoung])', 'r',6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'})
ylabel('Diff. P-HP Trial #')
[pval]=mynormandttest(xfplat([WTtouseyoung])-xfplathalf([WTtouseyoung]),xfplat([MUTtouseyoung])-xfplathalf([MUTtouseyoung]));
title([' p= ' num2str(pval)])
ylim([0 5000])
xlim([0 1.5])

% substraction from plateau and halfplateau 10-12mo 

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, xfplat([WTtouseold])'-xfplathalf([WTtouseold])', 'k',6,1);
paper_dotbarplot_modif(1, xfplat([MUTtouseold])'-xfplathalf([MUTtouseold])', 'r',6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'})
ylabel('Diff. P-HP Trial #')
[pval]=mynormandttest(xfplat([WTtouseold])-xfplathalf([WTtouseold]),xfplat([MUTtouseold])-xfplathalf([MUTtouseold]));
title([' p= ' num2str(pval)])
ylim([0 5000])
xlim([0 1.5])

%% Correlation of weight loss and performance
NormalizedWeights=(weights(:,4:end)*100)./weights(:,4);

% figure
% hold all
for ii=[1:8]
    for di=1:length(data(ii).expday)
        subplot(1,2,1)
        hold all
        plot(NormalizedWeights(ii,di),data(ii).dprime{1,di}(1), 'ko')
        title('Reinforced')
        subplot(1,2,2)
        hold all
        plot(NormalizedWeights(ii,di),data(ii).dprime{1,di}(2), 'ko')
        title('Probe')
    end
end


for ii=[1:8]
    for di=1:length(data(ii).expday)
        xw{ii}(di,:)=NormalizedWeights(ii,di);
        yrw{ii}(di,:)=data(ii).dprime{1,di}(1);
        ypw{ii}(di,:)=data(ii).dprime{1,di}(2);
    end
end

figure(1002)
subplot(1,2,1)
hold all
title('Reinforced')
for ii=MUTtouseyoung
plot(xw{ii},yrw{ii},'yo')
lsline
ylabel('d''')
xlabel('Weight %')
end
subplot(1,2,2)
hold all
title('Probe')
for ii=MUTtouseyoung
plot(xw{ii},ypw{ii},'yo')
lsline
end

%%

% figure
% subplot(2,1,1)
% hold all
% plot(xtrialsR(WTtouseyoung,:), dp_sli100R(WTtouseyoung, :), 'ko')
% plot(xtrialsR(MUTtouseyoung,:), dp_sli100R(MUTtouseyoung, :), 'ro')
% title('Reinforced')
% subplot(2,1,2)
% hold all
% plot(xtrialsP(WTtouseyoung,:), dp_slitouseP(WTtouseyoung, :), 'ko')
% plot(xtrialsP(MUTtouseyoung,:), dp_slitouseP(MUTtouseyoung, :), 'ro')
% abc=250;
% for ii=1:28
%     line([abc abc],[-1 5],'color', 'k')
%     abc=abc+250;
% end 
% xlabel('Mean trial #')
% ylabel('d''')
% title('Probe')

%% FIGURE 6: Tone exposure
% colorsmat2=[96, 211, 148; 238, 96, 85]./255; %[[155, 204, 85]./255; [228, 87, 46]./255; [179, 201, 174]./255; [239, 149, 137]./255];
colorsmat2=[96, 211, 148; 242, 143, 59]./255; 

figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0, NumTones(WTtouseyoung,1)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(0.5, NumTones(WTtouseyoung,2)', colorsmat2(2,:),6,1);
paper_dotbarplot_modif(1.5, NumTones(WTtouseyoung,1)', colorsmat2(1,:),6,1);
paper_dotbarplot_modif(2, NumTones(MUTtouseyoung,2)', colorsmat2(2,:),6,1);
set(gca, 'Xtick', [0.25 1.75], 'Xticklabel', {'Wt', 'APP^+'})
text([-0.4],[0.7], 'Target tone', 'Color', colorsmat2(1,:),'FontSize', 12, 'fontname','arial')
text([-0.4],[0.65], 'Foil tone', 'Color', colorsmat2(2,:),'FontSize', 12, 'fontname','arial')
ylabel('Probability')
ylim([0 0.75])
xlim([-0.5 2.5])
[h, p]=ranksum(NumTones(WTtouseyoung,1), NumTones(MUTtouseyoung,1));
[h, p]=ranksum(NumTones(MUTtouseyoung,1), NumTones(MUTtouseyoung,2));

%% Lucidity moments
lucidstart=NaN(32,50);
lucidsum=NaN(32,50);
lucidtotal=NaN(32,50);

i=0;
for ii=1:length(data)
    i=i+1;
%     clear lucidstart lucidsum lucidtotal
    daymaxPdprime(ii,:)=min(find(dprimePtoplot(:,ii)==max(dprimePtoplot(:,ii))));
    daymaxRdprime(ii,:)=min(find(dp_sli100R(:,ii)==max(dp_sli100R(:,ii))));
    for di=1:size(dp_lucid{ii},1)
        if ~isempty(find(dp_lucid{ii}(di,:)>max(dprimePtoplot(:,ii))))
            lucidstart(ii,di)=1;
            lucidsum(ii,di)=sum(dp_lucid{ii}(di,:)>max(dprimePtoplot(:,ii)));
            lucidtotal(ii,di)=sum(~isnan(dp_lucid{ii}(di,:)));
        else
            lucidstart(ii,di)=NaN;
            lucidsum(ii,di)=NaN;
            lucidtotal(ii,di)=NaN;
        end
    end
end

% Day of max probe and reinforced
figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, daymaxPdprime(WTtouseyoung,:)',colors(1,:),6,1);
paper_dotbarplot_modif(1.5, daymaxPdprime(MUTtouseyoung,:)',colors(2,:),6,1);
paper_dotbarplot_modif(1, daymaxRdprime(WTtouseyoung,:)', 'k',6,1);
paper_dotbarplot_modif(2, daymaxRdprime(MUTtouseyoung,:)', 'r',6,1);
[pval]=mynormandttest(daymaxPdprime(WTtouseyoung,:),daymaxPdprime(MUTtouseyoung,:));
[pval2]=mynormandttest(daymaxRdprime(WTtouseyoung,:),daymaxRdprime(MUTtouseyoung,:));
[pval3]=mynormandttest(daymaxPdprime(WTtouseyoung,:),daymaxRdprime(WTtouseyoung,:));
[pval4]=mynormandttest(daymaxPdprime(MUTtouseyoung,:),daymaxRdprime(MUTtouseyoung,:));
% ylim([-20 27])
xlim([0 2.5])
ylabel('Day of max d''')
line([0.5 1], [max(daymaxPdprime(MUTtouseyoung,:))+3 max(daymaxPdprime(MUTtouseyoung,:))+3], 'color', 'k');
line([1 2], [max(daymaxPdprime(MUTtouseyoung,:))+8 max(daymaxPdprime(MUTtouseyoung,:))+8], 'color', 'k');
line([0.5 1.5], [max(daymaxPdprime(MUTtouseyoung,:))+12 max(daymaxPdprime(MUTtouseyoung,:))+12], 'color', 'k');
line([1.5 2], [max(daymaxPdprime(MUTtouseyoung,:))+3 max(daymaxPdprime(MUTtouseyoung,:))+3], 'color', 'k');
txt1=(['p=', num2str(round(pval*100)/100)]);  txt2=(['p=', num2str(round(pval2*100)/100)]); txt3=(['p=', num2str(round(pval3*100)/100)]); txt4=(['p=', num2str(round(pval4*100)/100)]); 
text([0.5 0.5],[max(daymaxPdprime(MUTtouseyoung,:))+4 max(daymaxPdprime(MUTtouseyoung,:))+4],txt3, 'Color', 'k','FontSize', 12, 'fontname','arial');
text([1.5 1.5],[max(daymaxPdprime(MUTtouseyoung,:))+4 max(daymaxPdprime(MUTtouseyoung,:))+4],txt4, 'Color', 'k','FontSize', 12, 'fontname','arial');
text([0.5 0.5],[max(daymaxPdprime(MUTtouseyoung,:))+13 max(daymaxPdprime(MUTtouseyoung,:))+13],txt1, 'Color', 'k','FontSize', 12, 'fontname','arial');
text([1 1],[max(daymaxPdprime(MUTtouseyoung,:))+9 max(daymaxPdprime(MUTtouseyoung,:))+9],txt2, 'Color', 'k','FontSize', 12, 'fontname','arial');
set(gca, 'Xtick', [0.5 1 1.5 2], 'Xticklabel', {'Wt','APP^+','Wt','APP^+'});

% Number of "lucid" intrvals for 5 days starting when max probe d'
lucidprop=lucidsum./lucidtotal;
clear lucidpropaftdP
for ii=1:32
    lucidpropaftdP(ii,:)=nanmean(lucidprop(ii,daymaxPdprime(ii,:):daymaxPdprime(ii,:)+4));
end 

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, lucidpropaftdP(WTtouseyoung,:)','k',6);
paper_dotbarplot_modif(1, lucidpropaftdP(MUTtouseyoung,:)','r',6);
[pval]=mynormandttest(lucidpropaftdP(WTtouseyoung,:),lucidpropaftdP(MUTtouseyoung,:));
line([0.5 1], [max(lucidpropaftdP(MUTtouseyoung,:))+0.1 max(lucidpropaftdP(MUTtouseyoung,:))+0.1], 'color', 'k');
txt1=(['p=', num2str(round(pval*100)/100)]); 
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt','APP^+'});
xlim([0 1.5])
ylabel('Proportion of lucid / day')
text([0.5 0.5],[max(lucidpropaftdP(MUTtouseyoung,:))+0.15 max(lucidpropaftdP(MUTtouseyoung,:))+0.15],txt1, 'Color', 'k','FontSize', 12, 'fontname','arial');
ylim([0 1])

%% Lucidity moments with H and FA rates (dprime vs action rates 2 axis)
for di=[1:10]; %:10
    figure
    clear lucid_trans y
    i=0;
    for ii=[WTtouseyoung MUTtouseyoung]
        i=i+1; clear lucid_trans;
        subplot(4,4,i)
        hold all
        title([mouse{ii} ' d ' num2str(di)]); ylim([-1.5 4]);
        if size(dp_lucid{ii},1)>=di
            clear pd pdi
            yyaxis left
            hold all
            plot(dp_lucid{ii}(di,:), 'k-');
%             y=smooth([dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,:))))))) 0],5);
            line([1 400],[max(dprimePtouse{ii}) max(dprimePtouse{ii})],'color', 'm');
            xlim([0 300])
%             if sum(y>0)>1
%                 findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
%             [peaknum, peakloc]=findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
% 
%             end
            ylim([-1.5 6]); ylabel('d''');
            
            yyaxis right
            hold all
            plot(lucid_H{ii}(di,:), 'g-')
            plot(lucid_FA{ii}(di,:), 'r-')
            ylim([0 1.05]); ylabel('Action rate');
%             if ismember(ii, WTtouseyoung
            title([char(mouse(ii)) filesep num2str(di)])
            grid on
%                 title(['\color{Blue}' char(mouse(ii))])

        end
    end
end

%% Lucidity moments with H and FA rates (dprime vs criterion 2 axis)
for di=[1:10]; %:10
    figure
    clear lucid_trans y
    i=0;
    for ii=[WTtouseyoung MUTtouseyoung]
        i=i+1; clear lucid_trans;
        subplot(4,4,i)
        hold all
        title([mouse{ii} ' d ' num2str(di)]); ylim([-1.5 4]);
        if size(dp_lucid{ii},1)>=di
                        
            yyaxis right
            hold all
            patch([0 300 300 0],[-1.25 -1.25 1.25 1.25], [192,192,192]./255, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            plot(-1*c_lucid{ii}(di,:), 'r-')
            ylim([-3 3]); ylabel('-Criterion');xlim([0 300])
%             if ismember(ii, WTtouseyoung
            title([char(mouse(ii)) filesep num2str(di)])
            grid on
%                 title(['\color{Blue}' char(mouse(ii))])

            clear pd pdi
            yyaxis left
            hold all
            plot(dp_lucid{ii}(di,:), 'b-');
            %y=smooth([dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,:))))))) 0],5);
            line([1 400],[max(dprimePtouse{ii}) max(dprimePtouse{ii})],'color', 'm');
%             if sum(y>0)>1
%                 findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
%             [peaknum, peakloc]=findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
% 
%             end
            ylim([-1.5 6]); ylabel('d''');

        end
    end
end


%% Lucidity moments with H and FA rates
for di=1:10
    figure
    clear lucid_trans y
    i=0;
    for ii=[WTtouseyoung MUTtouseyoung]
        i=i+1; clear lucid_trans;
        subplot(2,4,i)
        hold all
        title([mouse{ii} ' d ' num2str(di)]); ylim([-1.5 4]);
        if size(dp_lucid{ii},1)>=di
            clear pd pdi
            yyaxis left
            plot(dp_lucid{ii}(di,:), 'k-');
            y=smooth([dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,:))))))) 0],5);
%             line([1 400],[max(dprimePtouse{ii}) max(dprimePtouse{ii})],'color', 'm');
%             if sum(y>0)>1
%                 findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
%             [peaknum, peakloc]=findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
% 
%             end
            ylim([-1.5 6]); ylabel('d''');
            
            yyaxis right
            plot(lucid_H{ii}(di,:), 'g-')
            plot(lucid_FA{ii}(di,:), 'r-')
            ylim([0 1.05]); ylabel('Action rate');
            title([char(mouse(ii)) filesep num2str(di)])
        end
    end
end

%%
for di=18
    figure
    clear lucid_trans y
    i=0;
    for ii=[WTtouseyoung MUTtouseyoung]
        i=i+1; clear lucid_trans;
        subplot(7,3,i)
        hold all
        title([mouse{ii} ' d ' num2str(di)]); ylim([-1.5 4]);
        if size(dp_lucid{ii},1)>=di
            clear pd pdi
            yyaxis left
            plot(dp_lucid{ii}(di,:), 'k-');
            y=smooth([dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,:))))))) 0],5);
            
            line([1 400],[max(dprimePtouse{ii}) max(dprimePtouse{ii})],'color', 'm');
            if sum(y>0)>1
                findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
                [peaknum, peakloc]=findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
                colors=[175, 252, 65; 0, 0, 0;58, 134, 255; 251, 86, 7;140, 138, 146;255, 190, 11]./255;
                
                behv_sav=[];
                for pi=1:length(peakloc)
                    plocH=nanmean(lucid_H{ii}(di,peakloc(pi)-2:peakloc(pi)+2));
                    plocFA=nanmean(lucid_FA{ii}(di,peakloc(pi)-2:peakloc(pi)+2));
                    
                    %TO OBTAIN MODES of lick behavior at different
                    %performance times
                    if  plocH>0.7 && plocFA<0.5 %if high hit and low fa showing high performance (GREAT INCREASE IN LEARNING)
                        behv_cat=1; 
                    elseif  plocH<0.5 && plocFA<0.5 %if both hit and fa decrease (DECREASE IN LICKING)
                        behv_cat=2;
                    elseif plocH>0.5 && plocFA>0.5 && plocH<=plocFA+(plocFA*0.1) %if both hit and fa increase (INCREASE IN LICKING)
                        behv_cat=3;
                    elseif  plocH>plocFA+(plocFA*0.1) && plocH>0.5 %if high hit and low fa showing high performance (MODERATE INCREASE IN LEARNING)
                        behv_cat=4;
                    elseif plocH>0.5 && plocFA>0.5 && abs(plocH-plocFA)>0.1 %slight increase in performance, only 0.1 difference above 0.5
                        behv_cat=6;
                    else
                        behv_cat=5;
                    end
                    behv_sav=[behv_sav behv_cat];
                    plot(peakloc(pi), peaknum(pi), 'o', 'color', colors(behv_cat,:), 'Markersize', 15,'LineWidth',3)
                    behv_cat_save{ii}(di).ilearn=sum(behv_sav==1);
                    behv_cat_save{ii}(di).dlick=sum(behv_sav==2);
                    behv_cat_save{ii}(di).ilick=sum(behv_sav==3);
                    behv_cat_save{ii}(di).mlearn=sum(behv_sav==4);
                    behv_cat_save{ii}(di).other=sum(behv_sav==5);
                    behv_cat_save{ii}(di).slearn=sum(behv_sav==6);
                end
                
            end
            ylim([-1.5 6]); ylabel('d''');
            yyaxis right
            plot(lucid_H{ii}(di,:), 'g-')
            plot(lucid_FA{ii}(di,:), 'r-')
            ylim([0 1.05]); ylabel('Action rate');
            title([char(mouse(ii)) filesep num2str(di)])
        end
    end
end

%% lucid intervals from plateau performance

clear trialnum ttrials
for ii=1:length(data)
    trialnum=cell2mat(data(ii).ntrials);
    for di=1:length(trialnum)
        ttrials{ii}(di)=nansum(trialnum(1:di));
    end
end

%%
% figure
i=0;
numlim=4;
dayaftplat=1;
clear sav_intH sav_intFA sav_breakp
for ii=[WTtouseyoung MUTtouseyoung]
    i=i+1;
%     subplot(7,3,i)
    [xvl,loc]=min(abs(ttrials{ii}-xfplat(ii)));
    loc=loc+dayaftplat;
%     locf=1:length(ttrials{ii});
    if length(ttrials{ii})>=loc+numlim
        locf=[loc:loc+numlim];
    else
        locf=[length(ttrials{ii})-numlim:length(ttrials{ii})];
    end
    
%     linespec={'-', '--', ':', '-.'};
%     hold all
    intH=[]; intFA=[]; clear breakpoint
    for si=1:length(locf)
        hnan=[]; fanan=[];
        hnan=find(~isnan(lucid_H{ii}(locf(si),:)));
        fanan=find(~isnan(lucid_FA{ii}(locf(si),:)));
        intH=[intH lucid_H{ii}(locf(si),hnan)];
        intFA=[intFA lucid_FA{ii}(locf(si),hnan)];
        if si>1
            breakpoint(si)=breakpoint(si-1)+length(lucid_H{ii}(locf(si),hnan));
        else
            breakpoint(si)=length(lucid_H{ii}(locf(si),hnan));
        end 
%         plot(lucid_H{ii}(locf(si),:), 'g', 'LineStyle', char(linespec(si)), 'linewidth', 2)
%         plot(lucid_FA{ii}(locf(si),:), 'r', 'LineStyle',char(linespec(si)), 'linewidth', 2)
%         plot(lucid_H{ii}(locf(si),:),lucid_FA{ii}(locf(si),:), 'o'); xlim([0 1]);
    end
    sav_intH{ii}=intH;
    sav_intFA{ii}=intFA;
    sav_breakp{ii}=breakpoint; 
%     ylim([0 1.05]); 
% %     ylabel('Action rate');
%     title([char(mouse(ii)) filesep num2str(locf)])
%     xlabel('H'); ylabel('FA');
%     legend({'4-ap'},'Location','northwest') %'1-ap',  '2-ap', '3-ap',
end

i=0;
figure
for ii=[MUTtouseyoung]
    i=i+1;
    subplot(length(MUTtouseyoung),1,i)
    hold all
    plot(sav_intH{ii}, 'g-')
    plot(sav_intFA{ii}, 'r-')
    for bip=1:length(sav_breakp{ii})
        line([sav_breakp{ii}(bip) sav_breakp{ii}(bip)],[0 1],'color', 'k', 'linestyle', ':');
    end
    ylabel('Probability')
    xlabel('Trials after plateau (5 days)')
    title([char(mouse(ii))])
    xlim([0 1600])
end

%% trial # vs std of fa and h
clear sav_st_intH sav_st_intFA winlength2
winlength=100;
i=0;
figure
for ii=[WTtouseyoung MUTtouseyoung]
    i=i+1;
    if i==11
        i=i+1;
    end
    
    wi=1;
    for ti=1:round((length(sav_intH{ii})-winlength)/winlength)
        winlength2{ii}(ti,:)=winlength*ti;
        sav_st_intH{ii}(ti,:)=nanstd(sav_intH{ii}(wi:winlength2{ii}(ti,:)));
        sav_st_intFA{ii}(ti,:)=nanstd(sav_intFA{ii}(wi:winlength2{ii}(ti,:)));
        wi=wi+winlength;
    end
    
    
    subplot(11,2,i)
    hold all
%     plot(winlength2{ii}, sav_st_intH{ii}, 'go')
    plot(winlength2{ii}, sav_st_intH{ii}, 'ro')
%     winlength2{ii}(isnan(sav_st_intFA{ii}))=[];
    sav_st_intH{ii}(isnan(sav_st_intH{ii}))=[];

    pfitFA(ii,:) = polyfit(winlength2{ii}, sav_st_intFA{ii},1);
    pfitH(ii,:) = polyfit(winlength2{ii}, sav_st_intH{ii},1);
    h=lsline; s=h.Color; h.Color='k'; t=h.LineWidth; h.LineWidth=2;
    ylabel('std p')
    xlabel('Trials after plateau (5 days)')
    title([char(mouse(ii))])
    xlim([0 8000])
    ylim([0 0.5])
end

%% Proportion of windows that fall higher than certain probability

clear prop_std_H prop_std_FA
psel=0.01;
for ii=[WTtouseyoung MUTtouseyoung]
    prop_std_H(ii)=sum(sav_st_intH{ii}>=psel)/length(sav_st_intH{ii});
    prop_std_FA(ii)=sum(sav_st_intFA{ii}>=psel)/length(sav_st_intFA{ii});
end 

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, prop_std_H(WTtouseyoung),'k',6,1);
paper_dotbarplot_modif(1, prop_std_H(MUTtouseyoung),'r',6,1);
[pval]=mynormandttest(prop_std_H(WTtouseyoung),prop_std_H(MUTtouseyoung));
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt','APP^+'});
xlim([0 1.5])
ylabel(['Proportion of std>' num2str(psel) ' pFAs'])
title(['p=' num2str(pval)])
% ylim([0 1])

%% slope of learning troughout for FA

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, pfitH(WTtouseyoung),'k',6,1);
paper_dotbarplot_modif(1, pfitH(MUTtouseyoung),'r',6,1);
[pval]=mynormandttest(pfitH(WTtouseyoung),pfitH(MUTtouseyoung));
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt','APP^+'});
xlim([0 1.5])
ylabel(['Learning slope for H'])
title(['p=' num2str(pval)])
% ylim([0 1])

%% Slope using per trial every 100 pfa ph

for ii=[WTtouseyoung MUTtouseyoung]
    
    xh=find(~isnan(pHsli100R(ii,:)));
    xfa=find(~isnan(pFAsli100R(ii,:)));
    sli100R_pfitH(ii,:) = polyfit(xtrialsR(ii,xh), pHsli100R(ii,xh),1);
    sli100R_pfitFA(ii,:) = polyfit(xtrialsR(ii,xfa), pFAsli100R(ii,xfa),1);
end

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, sli100R_pfitFA(WTtouseyoung),'k',6,1);
paper_dotbarplot_modif(1, sli100R_pfitFA(MUTtouseyoung),'r',6,1);
[pval]=mynormandttest(sli100R_pfitFA(WTtouseyoung),sli100R_pfitFA(MUTtouseyoung));
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt','APP^+'});
xlim([0 1.5])
ylabel(['Learning slope for H'])
title(['p=' num2str(pval)])

%% how much fa varies when h is high?
clear engaged_Diff


figure
i=0;
numlim=4; %how many days to consider in total calculation
dayaftplat=1; %how many days after plateau to start with
% ti_H=NaN(21,2000,1);
% ti_FA=NaN(21,2000,1);
for ii=[WTtouseyoung MUTtouseyoung]
    i=i+1;
    subplot(7,3,i)
    [xvl,loc]=min(abs(ttrials{ii}-xfplat(ii)));
    loc=loc+dayaftplat;
    
    if length(ttrials{ii})>=loc+numlim
        locf=[loc:loc+numlim];
    else
        locf=[length(ttrials{ii})-numlim:length(ttrials{ii})];
    end
    
%     plot(lucid_H{ii}(locf(si),:), lucid_FA{ii}(locf(si),:), 'o')

    hold all
    for si=1:length(locf)
        clear engagedHindx engagedH engagedFA 
        notengagedHindx=find(lucid_H{ii}(locf(si),:)>0.5 & lucid_H{ii}(locf(si),:)<0.95); %med engaged
        engagedHindx=find(lucid_H{ii}(locf(si),:)>=0.95); %engaged
%         engagedHindx=find(lucid_H{ii}(locf(si),:)<0.5); %not licking

        engagedH=lucid_H{ii}(locf(si),engagedHindx);
        notengagedH=lucid_H{ii}(locf(si),notengagedHindx);
        engagedFA=lucid_FA{ii}(locf(si),engagedHindx);
        notengagedFA=lucid_FA{ii}(locf(si),notengagedHindx);
        engaged_Diff(ii,si)=nanmean(engagedH-engagedFA);
        engaged_std(ii,si)=nanstd(engagedFA);
        notengaged_Diff(ii,si)=nanmean(notengagedH-notengagedFA);
        notengaged_std(ii,si)=nanstd(notengagedFA);
        maxmin=[max(engagedFA)-min(engagedFA)];
%         plot(engagedH, maxmin, 'o')
%         savH=lucid_H{ii}(locf(si),:); 
%         savFA=lucid_FA{ii}(locf(si),:);
%         savH(isnan(savH))=[];
%         savFA(isnan(savFA))=[];
        
%         ti_H(ii,:)=savH;
%         ti_FA(ii,:)=savFA;
        plot(lucid_H{ii}(locf(si),:), 'g-')
        plot(lucid_FA{ii}(locf(si),:), 'r-')
        plot(dp_lucid{ii}(locf(si),:), 'k-')
%         nanstd(dp_lucid{20}(23,:))/nanmean(dp_lucid{20}(23,:))
%         plot(lucid_H{ii}(locf(si),:), lucid_FA{ii}(locf(si),:), 'o')

%         plot(nanmean(engagedH-engagedFA), nanmean(notengagedH-notengagedFA), 'bo')
%         plot(engaged_std,engagedH-engagedFA, 'o', 'linewidth',2)
%         plot(engaged_Diff, engaged_std, 'o', 'linewidth',2)
        %         plot(lucid_H{ii}(locf(si),:),lucid_FA{ii}(locf(si),:), 'o'); xlim([0 1]);
    end
%     line([0 1], [0.1 0.1], 'linestyle', ':', 'color', 'k')
%     ylim([0 1]); xlim([0 1]);
    title([char(mouse(ii)) filesep num2str(locf)])
%     ylabel('std fa'); xlabel('diff h-fa');
%     xlabel('E h'); ylabel('E fa');
%     if ii==4
%         legend({'1-ap',  '2-ap', '3-ap','4-ap', '5-ap'},'Location','northwest') %'1-ap',  '2-ap', '3-ap',
%     end
end

%%

figure
hold all
plot(engaged_Diff(WTtouseyoung,1:11),engaged_std(WTtouseyoung,1:11), 'ko', 'linewidth',2)
plot(engaged_Diff(MUTtouseyoung,1:11),engaged_std(MUTtouseyoung,1:11), 'ro', 'linewidth',2)
line([0 1],[0.1 0.1], 'linestyle', ':', 'color', 'k')
xlim([0 1]); ylim([0 0.4]);
ylabel('std fa'); xlabel('diff h-fa');
title('Partially engaged in task')

figure
hold all
plot(notengaged_Diff(WTtouseyoung,1:11),notengaged_std(WTtouseyoung,1:11), 'ko', 'linewidth',2)
plot(notengaged_Diff(MUTtouseyoung,1:11),notengaged_std(MUTtouseyoung,1:11), 'ro', 'linewidth',2)
line([0 1],[0.1 0.1], 'linestyle', ':', 'color', 'k')
xlim([0 1]); ylim([0 0.4]);
ylabel('FA variability'); xlabel('Performance');



figure
hold all
plot(nanmean(engaged_std(WTtouseyoung,1:11),2),nanmean(engaged_Diff(WTtouseyoung,1:11),2), 'ko', 'linewidth',2)
plot(nanmean(engaged_std(MUTtouseyoung,1:11),2),nanmean(engaged_Diff(MUTtouseyoung,1:11),2), 'ro', 'linewidth',2)
line([0.1 0.1], [0 1], 'linestyle', ':', 'color', 'k')
ylim([0 1]); xlim([0 0.5]);
xlabel('std fa'); ylabel('diff h-fa');




%% Quantify behavior after plateau
clear dp_plateau
% endday=3;
maxnum=2;
for ii=[WTtouseyoung MUTtouseyoung]
    plateau_start=xfplat(ii);
    indxplateau=find(xtrialsR(ii, :)>plateau_start);
    yfit_dpval=yfit_dp_sli100R{ii}(2:end);
    platvalmax=sort(yfit_dpval(indxplateau:end));
    platvalmax(isnan(platvalmax))=[];
    if length(platvalmax)>=maxnum
        dp_plateau(ii,:)=nanmean(platvalmax(end-(maxnum-1):end));
    else
        dp_plateau(ii,:)=nanmean(platvalmax);
    end
end

clear dp_plateau
maxnum=5;
for ii=[WTtouseyoung MUTtouseyoung]
    plateau_start=xfplat(ii);
    indxplateau=find(xtrialsR(ii, :)>plateau_start);
    yfit_dpval=dp_sli100R(ii,indxplateau:end);
    platvalmax=sort(yfit_dpval);
    platvalmax(isnan(platvalmax))=[];
    if length(platvalmax)>=maxnum
        dp_plateau(ii,:)=nanmean(platvalmax(end-(maxnum-1):end));
    else
        dp_plateau(ii,:)=nanmean(platvalmax);
    end
end

clear dp_plateau
for ii=[WTtouseyoung MUTtouseyoung]
    plateau_start=xfplat(ii);
%     indxplateau=find(xtrialsR(ii, :)==plateau_start);
    yfit_dpval=dp_sli100R(ii,extraday);
    yfit_dpval=dp_sli100R(ii,:);
    dp_plateau(ii,:)=nanmean(yfit_dpval);
end

[a,b]=max(sum(xtrialsR>mean(xfplat)));
clear dp_plateau
for ii=[WTtouseyoung MUTtouseyoung]
    yfit_dpval=dp_sli100R(ii,b:end);
%     yfit_dpval=dp_sli100R(ii,:);
    dp_plateau(ii,:)=nanmean(yfit_dpval);
end


% fin

clear dp_plateau
% maxnum=5;
for ii=[WTtouseyoung MUTtouseyoung]
    plateau_start=5000;
    indxplateau=find(xtrialsR(ii, :)>plateau_start & xtrialsR(ii, :)<6000);
    yfit_dpval=dp_sli100R(ii,indxplateau:end);
%     platvalmax=sort(yfit_dpval);
%     platvalmax(isnan(platvalmax))=[];
%     if length(platvalmax)>=maxnum
        dp_plateau(ii,:)=nanmean(yfit_dpval);
%     else
%         dp_plateau(ii,:)=nanmean(platvalmax);
%     end
end

%Final figure REINFORCED
clear dp_plateau platvalmax
% endday=3;
maxnum=200;
for ii=1:16 %[WTtouseyoung MUTtouseyoung]
    plateau_start=xfplat(ii);
    indxplateau=find(xtrialsR(ii, :)>=plateau_start);
    yfit_dpval=yfit_dp_sli100R{ii}(2:end);
    platvalmax=sort(yfit_dpval(indxplateau:end));
    platvalmax(isnan(platvalmax))=[];
    if length(platvalmax)>=maxnum
        dp_plateau(ii,:)=nanmean(platvalmax(end-(maxnum-1):end));
    elseif length(platvalmax)==1
        dp_plateau(ii,:)=platvalmax;
    else
        dp_plateau(ii,:)=nanmean(platvalmax);
    end
end


figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateau(WTtouseyoung,:)', 'k',6,1);
paper_dotbarplot_modif(1, dp_plateau(MUTtouseyoung,:)', 'r',6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'})
ylabel('d''')
% ylim([-1 4])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateau(WTtouseyoung,:),dp_plateau(MUTtouseyoung,:));
title(['p=' num2str(pval)])

%old cohort
figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateau(WTtouseold,:)', 'k',6,1);
paper_dotbarplot_modif(1, dp_plateau(MUTtouseold,:)', 'r',6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'})
ylabel('d''')
ylim([-0.5 5])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateau(WTtouseold,:),dp_plateau(MUTtouseold,:));
title(['p=' num2str(pval)])

%Final figure PROBE
colors=[[187, 186, 183]; [255, 195, 185]]./255;
clear dp_plateauP platvalmax
% endday=3;
maxnum=200;
for ii=1:16; %[WTtouseyoung MUTtouseyoung]
    plateau_start=xfplat(ii);
    indxplateau=find(xtrialsP(ii, :)>=plateau_start); %xtrialsP
    yfit_dpval=yfit_dp_sli100P{ii}(2:end); % yfit_dp_sli100P
    platvalmax=sort(yfit_dpval(indxplateau:end));
    platvalmax(isnan(platvalmax))=[];
    if length(platvalmax)>=maxnum
        dp_plateauP(ii,:)=nanmean(platvalmax(end-(maxnum-1):end));
    elseif length(platvalmax)==1
        dp_plateauP(ii,:)=platvalmax;
    else
        dp_plateauP(ii,:)=nanmean(platvalmax);
    end
end


figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateauP(WTtouseyoung,:)', colors(1,:),6,1);
paper_dotbarplot_modif(1, dp_plateauP(MUTtouseyoung,:)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'})
ylabel('d''')
% ylim([0 5])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateauP(WTtouseyoung,:),dp_plateauP(MUTtouseyoung,:));
title(['p=' num2str(pval)])

%old cohort
figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateauP(WTtouseold,:)', colors(1,:),6,1);
paper_dotbarplot_modif(1, dp_plateauP(MUTtouseold,:)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'})
ylabel('d''')
ylim([-0.5 5])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateauP(WTtouseold,:),dp_plateauP(MUTtouseold,:));
title(['p=' num2str(pval)])

%COMPARE REINFORCED VS PROBE IN WT
figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateau(WTtouseyoung,:)', 'k',6,1);
paper_dotbarplot_modif(1, dp_plateauP(WTtouseyoung,:)', colors(1,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'R', 'P'})
ylabel('d''')
ylim([0 5])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateau(WTtouseyoung,:),dp_plateauP(WTtouseyoung,:));
title(['p=' num2str(pval)])


%COMPARE REINFORCED VS PROBE IN WT
figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateau(MUTtouseyoung,:)', 'r',6,1);
paper_dotbarplot_modif(1, dp_plateauP(MUTtouseyoung,:)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'R', 'P'})
ylabel('d''')
ylim([0 5])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateau(MUTtouseyoung,:),dp_plateauP(MUTtouseyoung,:));
title(['p=' num2str(pval)])


%COMPARE REINFORCED VS PROBE IN APP selective subpopulation bad
figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateau(indss,:)', 'r',6,1);
paper_dotbarplot_modif(1, dp_plateauP(indss,:)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'R', 'P'})
ylabel('d''')
ylim([0 5])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateau(indss,:),dp_plateauP(indss,:));
title(['p=' num2str(pval)])

%COMPARE REINFORCED VS PROBE IN APP selective subpopulation good
figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0.5, dp_plateau(new2MUTtouse,:)', 'r',6,1);
paper_dotbarplot_modif(1, dp_plateauP(new2MUTtouse,:)', colors(2,:),6,1);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'R', 'P'})
ylabel('d''')
ylim([0 5])
xlim([0 1.5])
[pval]=mynormandttest(dp_plateau(new2MUTtouse,:),dp_plateauP(new2MUTtouse,:));
title(['p=' num2str(pval)])


%% Distribution confidence internval
%tstatisticdprime interval
% intervalwt=[2.348878 3.294722];

intervalwt=[2.464358 3.303042]

indss=MUTtouseyoung(find(dp_plateau(MUTtouseyoung,:)<=intervalwt(1)));

new2MUTtouse=MUTtouseyoung(~ismember(MUTtouseyoung,indss));

[pval]=mynormandttest(dp_plateau(new2MUTtouse,:),dp_plateau(indss,:));


figure
hold all
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
paper_dotbarplot_modif(0, dp_plateau(WTtouseyoung,:)', 'k',6,1);
paper_dotbarplot_modif(0.5, dp_plateau(new2MUTtouse,:)', 'r',6,1);
paper_dotbarplot_modif(1, dp_plateau(indss,:)', 'r',6,1);
set(gca, 'Xtick', [0 0.5 1], 'Xticklabel', {'Wt','APP^+', 'APP^+diffdist'})
ylabel('d''')
ylim([0 5])
xlim([-0.5 1.5])
[pval]=mynormandttest(dp_plateau(indss,:),dp_plateau(new2MUTtouse,:));
title(['p=' num2str(pval)])

mouse(new2MUTtouse)
mouse(indss)

mouse(WTtouseyoung(find(dp_plateau(WTtouseyoung,:)>3)))
mouse(WTtouseyoung(find(dp_plateau(WTtouseyoung,:)<=3)))

%%
clear sum_il_Wt sum_ml_Wt sum_il_Mut sum_ml_Mut
i=0;
for ii=WTtouseyoung
    i=i+1;
    for di=[16:19]
        if size(behv_cat_save{1,ii},2)>=di
            if ~isempty(behv_cat_save{1,ii}(di).ilearn)
                il(di,:)=behv_cat_save{1,ii}(di).ilearn;
                ml(di,:)=behv_cat_save{1,ii}(di).mlearn;
            else
                il(di,:)=NaN;
                ml(di,:)=NaN;
            end
        end
    end
    sum_il_Wt(i,:)=nansum(il);
    sum_ml_Wt(i,:)=nansum(ml);
end


i=0;
for ii=MUTtouseyoung
    i=i+1;
    for di=[16:19]
        if size(behv_cat_save{1,ii},2)>=di
            if ~isempty(behv_cat_save{1,ii}(di).ilearn) && size(behv_cat_save{1,ii},2)>=di
                il(di,:)=behv_cat_save{1,ii}(di).ilearn;
                ml(di,:)=behv_cat_save{1,ii}(di).mlearn;
            else
                il(di,:)=NaN;
                ml(di,:)=NaN;
            end
        end
    end
    sum_il_Mut(i,:)=nansum(il);
    sum_ml_Mut(i,:)=nansum(ml);
end


figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, [sum_ml_Wt]','k',10);
paper_dotbarplot_modif(1, [sum_ml_Mut]','r',10);
xlim([0 1.5])
pval=mynormandttest([sum_ml_Wt], [sum_ml_Mut]);
ylabel('Peak learning events')
xlabel('Days 16-19')
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt','APP^+'});
title(['p= ' num2str(pval)])
ylim([0 50])

%%
dist_lucid_diff_wt=[];
dist_lucid_diff_mut=[];
clear lucid_trans lucid_diff
figure
i=0;
for di=1:4
    i=i+1;
    clear lucid_trans lucid_diff
    % figure
    for ii=[MUTtouseyoung WTtouseyoung]
        %     i=i+1;
        %     subplot(3,4,i)
        %     hold all
        if size(dp_lucid{ii},1)>=di
            clear pd pdi
            %         plot(dp_lucid{ii}(di,:), 'k-');
            y=smooth([dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,1:sum(~isnan(dp_lucid{ii}(di,:))))))) 0],5);
            if sum(y>0)>1
                [pdi, pd]=findpeaks(y,'MinPeakProminence',0.25,'MinPeakDistance',10);
                if ~isempty(pdi)
                                    lucid_trans(ii,i)=length(pd);
                    lucid_diff{ii}=pdi/max(dprimePtouse{ii}); %(max(dprimePtouse{ii})
                else
                    lucid_diff{ii}=NaN;
                end
            else
                                lucid_trans(ii,i)=NaN;
                lucid_diff{ii}=NaN;
            end
            %         xlim([1 sum(~isnan(dp_lucid{ii}(di,:)))+1]);
            %         ylim([-1 4])
        else
                        lucid_trans(ii,i)=NaN;
            lucid_diff{ii}=NaN;
        end
    end
    
    %   figure
    hold all
    for ii=WTtouseyoung
        plot(di, lucid_diff{ii}, 'ko')
        dist_lucid_diff_wt=[dist_lucid_diff_wt; lucid_diff{ii}];
    end
    for ii=MUTtouseyoung
        plot(di, lucid_diff{ii}, 'ro')
        dist_lucid_diff_mut=[dist_lucid_diff_mut; lucid_diff{ii}];
    end
    
%     [pd(ii,:),d]=ranksum()
    
end

%% Lucid closeness to max d' probe per animal


figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_beeswarm(0.5, dist_lucid_diff_wt','k',1);
paper_dotbarplot_beeswarm(1, dist_lucid_diff_mut','r',1);
xlim([0 1.5])
pval=mynormandttest(dist_lucid_diff_wt, dist_lucid_diff_mut);
ylabel('Max transitions trend')
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt','APP^+'});
title(['p= ' num2str(pval)])
ylim([-0.5 1.5])


% Figure 7: Amount of transitions per day
% figure
subplot(5,4,i)
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, lucid_trans(WTtouseyoung,:)','k',6);
paper_dotbarplot_modif(1, lucid_trans(MUTtouseyoung,:)','r',6);
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt','APP^+'});
xlim([0 1.5])
ylabel('# transitions / day ')
ylim([0 12])
title(['Day ' num2str(di)])

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, nanmean(lucid_trans(WTtouseyoung,:),2)','k',6);
paper_dotbarplot_modif(1, nanmean(lucid_trans(MUTtouseyoung,:),2)','r',6);
ylabel('# transitions / day (7-10)')


figure
hold all
for ii=WTtouseyoung
    plot(di, lucid_diff{ii}, 'ko')
end 
for ii=MUTtouseyoung
    plot(di, lucid_diff{ii}, 'ro')
end 
%% Lucidity variance analysis

for ii=1:32
    for di=1:size(dp_lucid{ii},1)
        var_dp_lucid{ii}(di,:)=nanstd(dp_lucid{ii}(di, :));
        mean_dp_lucid{ii}(di,:)=nanmean(dp_lucid{ii}(di, :));
    end
end

meanwt=NaN(35,32);
meanmut=NaN(35,32);

figure(15)
subplot(4,1,1)
hold all
title('WTyoung')
for ii=WTtouseyoung
  plot(var_dp_lucid{ii}./mean_dp_lucid{ii}, '-')  
end 
ylabel('cv')
subplot(4,1,2)
hold all
title('MUTyoung')
for ii=MUTtouseyoung
  plot(var_dp_lucid{ii}./mean_dp_lucid{ii}, '-')  
end 
ylabel('cv')
subplot(4,1,3)
hold all
for ii=WTtouseyoung
  plot(var_dp_lucid{ii}./mean_dp_lucid{ii}, 'k-')  
  meanwt(1:length(var_dp_lucid{ii}), ii)=var_dp_lucid{ii}./mean_dp_lucid{ii};
end 
for ii=MUTtouseyoung
  plot(var_dp_lucid{ii}./mean_dp_lucid{ii}, 'r-')  
  meanmut(1:length(var_dp_lucid{ii}), ii)=var_dp_lucid{ii}./mean_dp_lucid{ii};
end 
ylabel('cv')
subplot(4,1,4)
hold all
plot(nanmean(meanwt'), 'k-', 'linewidth', 3)
plot(nanmean(meanmut'), 'r-', 'linewidth', 3)
ylabel('Coefficient of variation')

%% FIGURE 3: Behavioral sliding d' variability
trialaprox1=NaN(32,50);
micetou=[WTtouseyoung MUTtouseyoung];

for ii=1:16
    clear ntrialsperday 
    ntrialsperday=cell2mat(data(ii).ntrials);
    for di=1:length(ntrialsperday)
        if length(ntrialsperday)>=d2
            trialaprox1(ii,di)=sum(ntrialsperday(1:di));
        else
            trialaprox1(ii,di)=NaN;
        end
    end
end

% Behavioral performance variability
[h, p1]=vartest2(nanmean(meanwt(1:4, :)), nanmean(meanmut(1:4, :)));
[h, p2]=vartest2(nanmean(meanwt(7:10, :)), nanmean(meanmut(7:10, :)));
[h, p3]=vartest2(nanmean(meanwt(16:19, :)), nanmean(meanmut(16:19, :)));

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, nanmean(meanwt(1:4, :)), 'k',6,1);
paper_dotbarplot_modif(1, nanmean(meanmut(1:4, :)), 'r',6,1);
paper_dotbarplot_modif(2, nanmean(meanwt(7:10, :)), 'k',6,1);
paper_dotbarplot_modif(2.5, nanmean(meanmut(7:10, :)), 'r',6,1);
paper_dotbarplot_modif(3.5, nanmean(meanwt(16:19, :)), 'k',6,1);
paper_dotbarplot_modif(4, nanmean(meanmut(16:19, :)), 'r',6,1);
ylim([-20 27])
xlim([0 4.5])
ylabel('Coeficient of variation')
set(gca, 'Xtick', [0.75 2.25 3.75], 'Xticklabel', {'~1-1000', '~2000-3000', '~5000-6000'});
text([3.5 3.5],[25 25], 'Wt', 'Color', 'k','FontSize', 12, 'fontname','arial')
text([3.5 3.5],[22 22], 'APP^+', 'Color', 'r','FontSize', 12, 'fontname','arial')
xlabel('Trials')
line([0.5 1], [max(nanmean(meanwt(1:4, :)))+5 max(nanmean(meanwt(1:4, :)))+5], 'color', 'k');
line([2 2.5],[max(nanmean(meanwt(7:10, :)))+5 max(nanmean(meanwt(7:10, :)))+5], 'color', 'k');
line([3.5 4],[max(nanmean(meanwt(16:19, :)))+5 max(nanmean(meanwt(16:19, :)))+5], 'color', 'k');
txt1=(['p=', num2str(round(p1*100)/100)]); txt2=(['p=', num2str(round(p2*100)/100)]); txt3=(['p=' num2str(p3)]);
text([0.5 0.5],[[max(nanmean(meanwt(1:4, :)))+5.7 max(nanmean(meanwt(1:4, :)))+5.7]],txt1, 'Color', 'k','FontSize', 12, 'fontname','arial');
text([2 2],[[max(nanmean(meanwt(7:10, :)))+5.7 max(nanmean(meanwt(7:10, :)))+5.7]],txt2, 'Color', 'k','FontSize', 12, 'fontname','arial');
text([3.5 3.5],[[max(nanmean(meanwt(16:19, :)))+5.7 max(nanmean(meanwt(16:19, :)))+5.7]],txt3, 'Color', 'k','FontSize', 12, 'fontname','arial');

%performance on those days of variability
[p,h]=ranksum(nanmean(dp_sli100R(WTtouseyoung,51:55),2),nanmean(dp_sli100R(WTtouseyoung,56:60),2));
[p,h]=ranksum(nanmean(dp_sli100R(MUTtouseyoung,51:55),2),nanmean(dp_sli100R(MUTtouseyoung,56:60),2));

dprangwt=(['d''=' num2str(nanmean(nanmean(dp_sli100R(WTtouseyoung,51:60),2)))]);
dprangmut=(['d''=' num2str(nanmean(nanmean(dp_sli100R(MUTtouseyoung,51:60),2)))]);

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(1, nanmean(meanwt(16:19, :)), 'k',6);
paper_dotbarplot_modif(1.5, nanmean(meanmut(16:19, :)), 'r',6);
ylim([0 6])
xlim([0.5 2])
ylabel('Coeficient of variation')
set(gca, 'Xtick', [1 1.5], 'Xticklabel', {'Wt', 'APP^+'});
text([0.55 0.55],[5.6 5.6], dprangwt, 'Color', 'k','FontSize', 12, 'fontname','arial')
text([0.55 0.55],[4.6 4.6], dprangmut, 'Color', 'r','FontSize', 12, 'fontname','arial')
xlabel('Trials ~5000-6000')
line([1 1.5],[max(nanmean(meanmut(16:19, :)))+0.5 max(nanmean(meanmut(16:19, :)))+0.5], 'color', 'k');
txt3=(['p=' num2str(p3)]);
text([1 1],[[max(nanmean(meanmut(16:19, :)))+1 max(nanmean(meanmut(16:19, :)))+1]],txt3, 'Color', 'k','FontSize', 12, 'fontname','arial');

%% coef variation after plateau
i=0;
numlim=4; %how many days to consider in total calculation
dayaftplat=1; %how many days after plateau to start with
for ii=[WTtouseyoung]
    i=i+1;
    [xvl,loc]=min(abs(ttrials{ii}-xfplat(ii)));
    loc=loc+dayaftplat;
    
    if length(ttrials{ii})>=loc+numlim
        locf=[loc:loc+numlim];
    else
        locf=[length(ttrials{ii})-numlim:length(ttrials{ii})];
    end
    
    var_plat_wt(i,:)=nanmean(meanwt(locf,ii));
end 


i=0;
numlim=4; %how many days to consider in total calculation
dayaftplat=1; %how many days after plateau to start with
for ii=[MUTtouseyoung]
    i=i+1;
    [xvl,loc]=min(abs(ttrials{ii}-xfplat(ii)));
    loc=loc+dayaftplat;
    
    if length(ttrials{ii})>=loc+numlim
        locf=[loc:loc+numlim];
    else
        locf=[length(ttrials{ii})-numlim:length(ttrials{ii})];
    end
    
    var_plat_mut(i,:)=nanmean(meanmut(locf,ii));
end 

figure
set(gca,'TickDir','out','FontSize',12,'fontname','arial');
hold all
paper_dotbarplot_modif(0.5, var_plat_wt', 'k',6,1);
paper_dotbarplot_modif(1, var_plat_mut', 'r',6,1);
pval=mynormandttest(var_plat_wt, var_plat_mut);
title(['p=' num2str(pval)])
% ylim([0 6])
xlim([0 1.5])
ylabel('Coeficient of variation')
set(gca, 'Xtick', [0.5 1], 'Xticklabel', {'Wt', 'APP^+'});

%% Distribution of sliding window dprimes 

bins=[-2:0.25:5];

for ii=1:32
    for di=1:size(dp_lucid{ii},1)
        for bi=1:length(bins)-1
            dp_binned{ii}(di,bi)=sum(dp_lucid{ii}(di,:)>bins(bi) & dp_lucid{ii}(di,:)<=bins(bi+1))/sum(~isnan(dp_lucid{ii}(di,:)));
            binac(bi,:)=bins(bi)+(bins(bi+1)-bins(bi))/2;
        end
%         daymaxPdprime=min(find(dprimePtoplot(:,ii)==max(dprimePtoplot(:,ii))));
%         figure(ii)
%         set(gcf, 'Position', get(0, 'Screensize'));
%         subplot(6,6,di)
%         hold all
%         set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
%         plot(dp_binned{ii}(di,:), 'k-', 'linewidth', 2)
%         set(gca, 'Xtick', [1:4:length(bins)], 'Xticklabel', [bins(1:4:end)])
%         ylim([0 1])
%         xlabel('d''')
%         ylabel('Proportion')
%         title([char(mouse(ii)) filesep num2str(di)])
%         line([find(bins==1) find(bins==1)],[0 1], 'color', 'k', 'linestyle', ':');
%         if di>=daymaxPdprime
%            indnear=interp1(bins,1:length(bins),max(dprimePtoplot(:,ii)),'nearest');
%            line([indnear indnear],[0 1], 'color', 'b', 'linestyle', '-');
%         end 
    end
end

%plot when distribution reaches above 1 clearly
wtdist=[14, 11, 13, 6, 12, 12, 13, 13, 12, 14];
mutdist=[19, 17, 8, 5, 17, 13, 15]; %4 mice did never reach that distribution

[h, p]=ranksum(wtdist, mutdist);

% width
% dy=9;
i=0;
for dy=1:20
    clear widthdp
    for ii=[WTtouseyoung MUTtouseyoung];
        i=i+1; 
        if size(dp_binned{ii}, 1)>=dy
            widthdp(ii, :)=fwhm([binac'], dp_binned{ii}(dy,:));
        else
            widthdp(ii, :)=NaN;
        end 
        %     figure(10002)
        %     subplot(6,6,i)
        %     hold all
        %     plot(dp_binned{ii}(dy,:))
        %     title([char(mouse(ii)) filesep num2str(dy) filesep num2str(widthdp(ii,:))])
    end
    figure(190)
    subplot(4,5,dy)
    hold all
    [out]=paper_dotbarplot_modif(1, widthdp(WTtouseyoung,:)', 'k')
    [out]=paper_dotbarplot_modif(1+1, widthdp(MUTtouseyoung,:)', 'r')
    [p, h]=ranksum(widthdp(WTtouseyoung,:), widthdp(MUTtouseyoung,:));
    title(['day ' num2str(dy) filesep ' p=' num2str(p)])
    ylim([0 2.5])
    ylabel('fwhm')
    set(gca, 'Xtick', [1 2], 'Xticklabel', {'wt', 'app+'})
end

%%

%young
xRWT=nanmean(xsaveR(:, WTtouseyoung),2);
xRMUT=nanmean(xsaveR(:, MUTtouseyoung2),2);
xPWT=nanmean(xsaveP(:, WTtouseyoung),2);
xPMUT=nanmean(xsaveP(:, MUTtouseyoung2),2);

[param ,statxRWT]=sigm_fit([1:xi], xRWT');      
[param ,statxRMUT]=sigm_fit([1:xi], xRMUT');      
[param ,statxPWT]=sigm_fit([1:xi-2], xPWT(3:xi)');      
[param ,statxPMUT]=sigm_fit([1:xi-2], xPMUT(3:xi)');  

%old
xRWTo=nanmean(xsaveR(:, WTtouseold),2);
xRMUTo=nanmean(xsaveR(:, MUTtouseold),2);
xPWTo=nanmean(xsaveP(:, WTtouseold),2);
xPMUTo=nanmean(xsaveP(:, MUTtouseold),2);

[param ,statxRWTo]=sigm_fit([1:xi], xRWTo');      
[param ,statxRMUTo]=sigm_fit([1:xi], xRMUTo');      
[param ,statxPWTo]=sigm_fit([1:xi-2], xPWTo(3:xi)');      
[param ,statxPMUTo]=sigm_fit([1:xi-2], xPMUTo(3:xi)');  

%sem
xr=xsaveR;
xp=xsaveP;

SEMRY_WT=nanstd(xr(:,WTtouseyoung)')/sqrt(length(WTtouseyoung));
SEMRY_MUT=nanstd(xr(:,MUTtouseyoung)')/sqrt(length(MUTtouseyoung));
SEMPY_WT=nanstd(xp(:,WTtouseyoung)')/sqrt(length(WTtouseyoung));
SEMPY_MUT=nanstd(xp(:,MUTtouseyoung)')/sqrt(length(MUTtouseyoung));

SEMRO_WT=nanstd(xr(:,WTtouseold)')/sqrt(length(WTtouseold));
SEMRO_MUT=nanstd(xr(:,MUTtouseold)')/sqrt(length(MUTtouseold));
SEMPO_WT=nanstd(xp(:,WTtouseold)')/sqrt(length(WTtouseold));
SEMPO_MUT=nanstd(xp(:,MUTtouseold)')/sqrt(length(MUTtouseold));
    

%% D PRIME FINAL FIGS PER DAY
figure
subplot(2,3,1)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
plot([1:xi], nanmean(xsaveR(:, WTtouseyoung),2), 'ko', 'Markerfacecolor', 'k')
plot([1:xi], nanmean(xsaveR(:, MUTtouseyoung2),2), 'ro', 'Markerfacecolor', 'r')
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:end, WTtouseyoung),2)], 'o','color', colors(1, 1:3), 'Markerfacecolor', colors(1, 1:3))
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:end, MUTtouseyoung2),2)], 'o', 'color', colors(2, 1:3), 'Markerfacecolor', colors(2, 1:3))
plot([3:xi], statxPWT.ypred, '-', 'color', colors(1, 1:3), 'linewidth', 2);
plot([3:xi], statxPMUT.ypred, '-', 'color', colors(2, 1:3), 'linewidth', 2);
plot([1:xi], statxRWT.ypred, 'k-', 'linewidth', 2);
plot([1:xi], statxRMUT.ypred, 'r-', 'linewidth', 2);
line([1 xi+1], [1 1], 'color', 'k', 'linestyle', ':')
xlim([0 xi+1])
ylim([-1 4])
ylabel('d''')
xlabel('Training day')
legend(['Reinf Wt (',num2str(length(WTtouseyoung)) ')'], ['Reinf APP+ (',num2str(length(MUTtouseyoung2)) ')'],'Probe Wt', 'Probe APP+ ', 'location', 'Northwest')
title('Young')

subplot(2,3,4)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
plot([1:xi], nanmean(xsaveR(:, WTtouseold),2), 'ko', 'Markerfacecolor', 'k')
plot([1:xi], nanmean(xsaveR(:, MUTtouseold),2), 'ro', 'Markerfacecolor', 'r')
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:end, WTtouseold),2)], 'o','color', colors(1, 1:3), 'Markerfacecolor', colors(1, 1:3))
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:end, MUTtouseold),2)], 'o', 'color', colors(2, 1:3), 'Markerfacecolor', colors(2, 1:3))
plot([3:xi], statxPWTo.ypred, '-', 'color', colors(1, 1:3), 'linewidth', 2);
plot([3:xi], statxPMUTo.ypred, '-', 'color', colors(2, 1:3), 'linewidth', 2);
plot([1:xi], statxRWTo.ypred, 'k-', 'linewidth', 2);
plot([1:xi], statxRMUTo.ypred, 'r-', 'linewidth', 2);
line([1 xi+1], [1 1], 'color', 'k', 'linestyle', ':')
xlim([0 xi+1])
ylim([-1 4])
ylabel('d''')
xlabel('Training day')
legend(['Reinf Wt (',num2str(length(WTtouseold)) ')'], ['Reinf APP+ (',num2str(length(MUTtouseold)) ')'],'Probe Wt', 'Probe APP+ ', 'location', 'Northwest')
title('Old')

subplot(2,3,2)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
errorbar(nanmean(xsaveR(:, WTtouseyoung),2), SEMRY_WT, 'k', 'LineStyle','none');
errorbar(nanmean(xsaveR(:, MUTtouseyoung2),2), SEMRY_MUT, 'r', 'LineStyle','none');
plot([1:xi], nanmean(xsaveR(:, WTtouseyoung),2), 'ko', 'Markerfacecolor', 'k')
plot([1:xi], nanmean(xsaveR(:, MUTtouseyoung2),2), 'ro', 'Markerfacecolor', 'r')
plot([1:xi], statxRWT.ypred, 'k-', 'linewidth', 2);
plot([1:xi], statxRMUT.ypred, 'r-', 'linewidth', 2);
line([1 xi+1], [1 1], 'color', 'k', 'linestyle', ':')
xlim([0 xi+1])
ylim([-1 4])
ylabel('d''')
xlabel('Training day')
title('Reinforced')

subplot(2,3,3)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
errorbar([NaN;NaN; nanmean(xsaveP(3:xi, WTtouseyoung),2)], [NaN, NaN, SEMPY_WT(3:xi)], 'color', colors(1, 1:3), 'LineStyle','none');
errorbar([NaN;NaN; nanmean(xsaveP(3:xi, MUTtouseyoung2),2)], [NaN, NaN, SEMPY_MUT(3:xi)], 'color', colors(2, 1:3), 'LineStyle','none');
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:xi, WTtouseyoung),2)], 'o','color', colors(1, 1:3), 'Markerfacecolor', colors(1, 1:3))
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:xi, MUTtouseyoung2),2)], 'o', 'color', colors(2, 1:3), 'Markerfacecolor', colors(2, 1:3))
plot([3:xi], statxPWT.ypred, '-', 'color', colors(1, 1:3), 'linewidth', 2);
plot([3:xi], statxPMUT.ypred, '-', 'color', colors(2, 1:3), 'linewidth', 2);
line([1 xi+1], [1 1], 'color', 'k', 'linestyle', ':')
xlim([0 xi+1])
ylim([-1 4])
ylabel('d''')
xlabel('Training day')
title('Probe')

subplot(2,3,5)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
errorbar(nanmean(xsaveR(:, WTtouseold),2), SEMRO_WT, 'k', 'LineStyle','none');
errorbar(nanmean(xsaveR(:, MUTtouseold),2), SEMRO_MUT, 'r', 'LineStyle','none');
plot([1:xi], nanmean(xsaveR(:, WTtouseold),2), 'ko', 'Markerfacecolor', 'k')
plot([1:xi], nanmean(xsaveR(:, MUTtouseold),2), 'ro', 'Markerfacecolor', 'r')
plot([1:xi], statxRWTo.ypred, 'k-', 'linewidth', 2);
plot([1:xi], statxRMUTo.ypred, 'r-', 'linewidth', 2);
line([1 xi+1], [1 1], 'color', 'k', 'linestyle', ':')
xlim([0 xi+1])
ylim([-1 4])
ylabel('d''')
xlabel('Training day')
title('Reinforced')

subplot(2,3,6)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
errorbar([NaN;NaN; nanmean(xsaveP(3:xi, WTtouseold),2)], [NaN, NaN, SEMPO_WT(3:xi)], 'color', colors(1, 1:3), 'LineStyle','none');
errorbar([NaN;NaN; nanmean(xsaveP(3:xi, MUTtouseold),2)], [NaN, NaN, SEMPO_MUT(3:xi)], 'color', colors(2, 1:3), 'LineStyle','none');
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:xi, WTtouseold),2)], 'o','color', colors(1, 1:3), 'Markerfacecolor', colors(1, 1:3))
plot([1:xi], [NaN;NaN; nanmean(xsaveP(3:xi, MUTtouseold),2)], 'o', 'color', colors(2, 1:3), 'Markerfacecolor', colors(2, 1:3))
plot([3:xi], statxPWTo.ypred, '-', 'color', colors(1, 1:3), 'linewidth', 2);
plot([3:xi], statxPMUTo.ypred, '-', 'color', colors(2, 1:3), 'linewidth', 2);
line([1 xi+1], [1 1], 'color', 'k', 'linestyle', ':')
xlim([0 xi+1])
ylim([-1 4])
ylabel('d''')
xlabel('Training day')
title('Probe')

%% analysis dprime
tou=xsaveR;
anovaWT=WTtouseyoung;
anovaMUT=MUTtouseyoung;
clear anov_x anov_y anov_g anov_d
anov_x=[tou(1:11, anovaWT), tou(1:11, anovaMUT)];
anov_y=reshape(anov_x,1,[]); 
anov_g=[ones(1,length(anovaWT)*size(anov_x,1)), ones(1,length(anovaMUT)*size(anov_x,1))*2];
anov_d=repmat([1:size(anov_x,1)], 1, length(anov_y)/size(anov_x,1));
subj=[repmat([anovaWT], 1, size(anov_x,1)) repmat([anovaMUT], 1, size(anov_x,1))];

[p,tbl,stats]=anovan(anov_y, {anov_g, anov_d},'varnames',{'group','day'}, 'model','interaction');


stats = rm_anova2(Y,S,F1,F2,FACTNAMES);

stats = rm_anova2(anov_y',subj',anov_g',anov_d',{'genotype', 'day'});


%% hit and fa rate for days
daystouse=[1 2 5 10 15 21];
wtxpoint=[0:2:length(daystouse)*2]+0.7;
mutxpoint=wtxpoint+0.6;

figure
subplot(2,2,1)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
for ii=1:length(daystouse)
paper_dotbarplot_modif(wtxpoint(ii), hitRtoplot(daystouse(ii), WTtouseyoung)*100, 'k');
paper_dotbarplot_modif(mutxpoint(ii), hitRtoplot(daystouse(ii), MUTtouseyoung)*100, 'r');
line([wtxpoint(ii)-0.2 mutxpoint(ii)+0.2],[110 110], 'color', 'k');
clear pval txt
pval=mynormandttest(hitRtoplot(daystouse(ii), WTtouseyoung)*100, hitRtoplot(daystouse(ii), MUTtouseyoung)*100);
txt=(['p=' num2str(round(pval*100)/100)]);
text([wtxpoint(ii)-0.2],[115], txt);
end 
title('Young')
ylabel('Hit rate (%)')
xlabel('Training day')
set(gca, 'Xtick', round(wtxpoint(1:length(daystouse))), 'Xticklabel', daystouse);
xlim([0 max(mutxpoint(1:length(daystouse)))+0.7])
ylim([0 120])

subplot(2,2,2)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
for ii=1:length(daystouse)
paper_dotbarplot_modif(wtxpoint(ii), faRtoplot(daystouse(ii), WTtouseyoung)*100, 'k');
paper_dotbarplot_modif(mutxpoint(ii), faRtoplot(daystouse(ii), MUTtouseyoung)*100, 'r');
line([wtxpoint(ii)-0.2 mutxpoint(ii)+0.2],[110 110], 'color', 'k');
clear pval txt
pval=mynormandttest(faRtoplot(daystouse(ii), WTtouseyoung)*100, faRtoplot(daystouse(ii), MUTtouseyoung)*100);
txt=(['p=' num2str(round(pval*100)/100)]);
text([wtxpoint(ii)-0.2],[115], txt);
end 
title('Young')
ylabel('False alarm rate (%)')
xlabel('Training day')
set(gca, 'Xtick', round(wtxpoint(1:length(daystouse))), 'Xticklabel', daystouse);
xlim([0 max(mutxpoint(1:length(daystouse)))+0.7])
ylim([0 120])

subplot(2,2,3)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
for ii=1:length(daystouse)
paper_dotbarplot_modif(wtxpoint(ii), hitRtoplot(daystouse(ii), WTtouseold)*100, 'k');
paper_dotbarplot_modif(mutxpoint(ii), hitRtoplot(daystouse(ii), MUTtouseold)*100, 'r');
line([wtxpoint(ii)-0.2 mutxpoint(ii)+0.2],[110 110], 'color', 'k');
clear pval txt
pval=mynormandttest(hitRtoplot(daystouse(ii), WTtouseold)*100, hitRtoplot(daystouse(ii), MUTtouseold)*100);
txt=(['p=' num2str(round(pval*100)/100)]);
text([wtxpoint(ii)-0.2],[115], txt);
end 
title('Old')
ylabel('Hit rate (%)')
xlabel('Training day')
set(gca, 'Xtick', round(wtxpoint(1:length(daystouse))), 'Xticklabel', daystouse);
xlim([0 max(mutxpoint(1:length(daystouse)))+0.7])
ylim([0 120])

subplot(2,2,4)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
for ii=1:length(daystouse)
paper_dotbarplot_modif(wtxpoint(ii), faRtoplot(daystouse(ii), WTtouseold)*100, 'k');
paper_dotbarplot_modif(mutxpoint(ii), faRtoplot(daystouse(ii), MUTtouseold)*100, 'r');
line([wtxpoint(ii)-0.2 mutxpoint(ii)+0.2],[110 110], 'color', 'k');
clear pval txt
pval=mynormandttest(faRtoplot(daystouse(ii), WTtouseold)*100, faRtoplot(daystouse(ii), MUTtouseold)*100);
txt=(['p=' num2str(round(pval*100)/100)]);
text([wtxpoint(ii)-0.2],[115], txt);
end 
title('Old')
ylabel('False alarm rate (%)')
xlabel('Training day')
set(gca, 'Xtick', round(wtxpoint(1:length(daystouse))), 'Xticklabel', daystouse);
xlim([0 max(mutxpoint(1:length(daystouse)))+0.7])
ylim([0 120])

%% First lick latency per day
daystouse=[1:15]; %[1 2 5 10 15 21];
wtxpoint=[0:2:length(daystouse)*2]+0.7;
mutxpoint=wtxpoint+0.6;

% RTconfidenttoplot

figure
subplot(1,2,1)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
for ii=1:length(daystouse)
paper_dotbarplot_modif(wtxpoint(ii), RTconfidenttoplot(daystouse(ii), WTtouseyoung), 'b',6,1);
paper_dotbarplot_modif(mutxpoint(ii), RTconfidenttoplot(daystouse(ii), MUTtouseyoung), 'k',6,1);
line([wtxpoint(ii)-0.2 mutxpoint(ii)+0.2],[1 1], 'color', 'k');
clear pval txt
pval=mynormandttest(RTconfidenttoplot(daystouse(ii), WTtouseyoung), RTconfidenttoplot(daystouse(ii), MUTtouseyoung));
txt=(['p=' num2str(round(pval*100)/100)]);
text([wtxpoint(ii)-0.2],[1.05], txt);
end 
title('Young')
ylabel('Proportion of confident 1st lick (Target)')
xlabel('Training day')
set(gca, 'Xtick', round(wtxpoint(1:length(daystouse))), 'Xticklabel', daystouse);
xlim([0 max(mutxpoint(1:length(daystouse)))+0.7])
ylim([0 1])

subplot(1,2,2)
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
hold all
for ii=1:length(daystouse)
paper_dotbarplot_modif(wtxpoint(ii), RFconfidenttoplot(daystouse(ii), WTtouseyoung), 'b',6,1);
paper_dotbarplot_modif(mutxpoint(ii), RFconfidenttoplot(daystouse(ii), MUTtouseyoung), 'k',6,1);
line([wtxpoint(ii)-0.2 mutxpoint(ii)+0.2],[1 1], 'color', 'k');
clear pval txt
pval=mynormandttest(RFconfidenttoplot(daystouse(ii), WTtouseyoung), RFconfidenttoplot(daystouse(ii), MUTtouseyoung));
txt=(['p=' num2str(round(pval*100)/100)]);
text([wtxpoint(ii)-0.2],[1.05], txt);
end 
title('Young')
ylabel('Proportion of confident 1st lick (Foil)')
xlabel('Training day')
set(gca, 'Xtick', round(wtxpoint(1:length(daystouse))), 'Xticklabel', daystouse);
xlim([0 max(mutxpoint(1:length(daystouse)))+0.7])
ylim([0 1])




%% STatistics

[h, pwt]=kstest(xsaveR(:, WTtouseyoung));
[h, pmut]=kstest(xsaveR(:, MUTtouseyoung));

[p,tbl,stats] = friedman(xsaveR(1:10, WTtouseyoung),'on')

mouse(WTtouseyoung)

genotypename=[repmat({'wt'}, 1, length(WTtouseyoung)) repmat({'mut'}, 1, length(MUTtouseyoung))];



t1=table(genotypename', [xsaveR(1, WTtouseyoung)'; xsaveR(1, MUTtouseyoung)'], [xsaveR(2, WTtouseyoung)'; xsaveR(2, MUTtouseyoung)'], ...
    [xsaveR(3, WTtouseyoung)'; xsaveR(3, MUTtouseyoung)'], [xsaveR(4, WTtouseyoung)'; xsaveR(4, MUTtouseyoung)'], [xsaveR(5, WTtouseyoung)'; xsaveR(5, MUTtouseyoung)'],...
    [xsaveR(6, WTtouseyoung)'; xsaveR(6, MUTtouseyoung)'], [xsaveR(7, WTtouseyoung)'; xsaveR(7, MUTtouseyoung)'], ...
    [xsaveR(8, WTtouseyoung)'; xsaveR(8, MUTtouseyoung)'], [xsaveR(9, WTtouseyoung)'; xsaveR(9, MUTtouseyoung)'], [xsaveR(10, WTtouseyoung)'; xsaveR(10, MUTtouseyoung)'],...
    [xsaveR(11, WTtouseyoung)'; xsaveR(11, MUTtouseyoung)'],'VariableNames',{'genotype','meas1','meas2','meas3','meas4','meas5','meas6','meas7','meas8','meas9','meas10', 'meas11'});

t2=table(genotypename', [xsaveR(12, WTtouseyoung)'; xsaveR(12, MUTtouseyoung)'], ...
    [xsaveR(13, WTtouseyoung)'; xsaveR(13, MUTtouseyoung)'], [xsaveR(14, WTtouseyoung)'; xsaveR(14, MUTtouseyoung)'], [xsaveR(15, WTtouseyoung)'; xsaveR(15, MUTtouseyoung)'],...
    [xsaveR(16, WTtouseyoung)'; xsaveR(16, MUTtouseyoung)'], [xsaveR(17, WTtouseyoung)'; xsaveR(17, MUTtouseyoung)'], ...
    [xsaveR(18, WTtouseyoung)'; xsaveR(18, MUTtouseyoung)'], [xsaveR(19, WTtouseyoung)'; xsaveR(19, MUTtouseyoung)'], [xsaveR(20, WTtouseyoung)'; xsaveR(20, MUTtouseyoung)'],...
    [xsaveR(21, WTtouseyoung)'; xsaveR(21, MUTtouseyoung)'], 'VariableNames',{'genotype','meas12','meas13','meas14','meas15','meas16','meas17','meas18','meas19','meas20','meas21'});

t3=table(genotypename', [xsaveR(15, WTtouseyoung)'; xsaveR(15, MUTtouseyoung)'],...
    [xsaveR(16, WTtouseyoung)'; xsaveR(16, MUTtouseyoung)'], [xsaveR(17, WTtouseyoung)'; xsaveR(17, MUTtouseyoung)'], ...
    [xsaveR(18, WTtouseyoung)'; xsaveR(18, MUTtouseyoung)'], [xsaveR(19, WTtouseyoung)'; xsaveR(19, MUTtouseyoung)'], [xsaveR(20, WTtouseyoung)'; xsaveR(20, MUTtouseyoung)'],...
    [xsaveR(21, WTtouseyoung)'; xsaveR(21, MUTtouseyoung)'], 'VariableNames',{'genotype','meas15','meas16','meas17','meas18','meas19','meas20','meas21'});
clear MEAS rm ranovatbl
MEAS = table([12:21]','VariableNames',{'Measurements'});
rm = fitrm(t2,'meas12-meas21~genotype','WithinDesign',MEAS);
ranovatbl = ranova(rm)


mtx=[xsaveR(1:10, WTtouseyoung)'; xsaveR(1:10, MUTtouseyoung)'];
%%








%% FINAL Figures

colors=[[187, 186, 183]; [255, 195, 185]]./255;
endday=18;
figure
subplot(1,2,1)
hold all
errorbar(nanmean(dprimeRtoplot(1:endday,WTtouseyoung),2), SEMRY_WT(1:endday), 'k', 'linewidth', 2)
errorbar(nanmean(dprimeRtoplot(1:endday,MUTtouseyoung),2), SEMRY_MUT(1:endday), 'color', [229, 61, 45]./255, 'linewidth', 2)
errorbar(nanmean(dprimePtoplot(1:endday,WTtouseyoung),2), SEMPY_WT(1:endday), 'k:', 'linewidth', 2)
errorbar(nanmean(dprimePtoplot(1:endday,MUTtouseyoung),2), SEMPY_WT(1:endday), ':','color', [229, 61, 45]./255, 'linewidth', 2)
ylabel('d''')
xlabel('Day')
ylim([-1 6])
xlim([0 19])
legend(['Wt ',num2str(length(WTtouseyoung))], ['APP+ ',num2str(length(MUTtouseyoung))])
title('Young')

subplot(1,2,2)
hold all
errorbar(nanmean(dprimeRtoplot(1:endday,WTtouseold),2), SEMRO_WT(1:endday), 'k', 'linewidth', 2)
errorbar(nanmean(dprimeRtoplot(1:endday,MUTtouseold),2), SEMRO_MUT(1:endday), 'color', [229, 61, 45]./255, 'linewidth', 2)
errorbar(nanmean(dprimePtoplot(1:endday,WTtouseold),2), SEMPO_WT(1:endday), 'k:', 'linewidth', 2)
errorbar(nanmean(dprimePtoplot(1:endday,MUTtouseold),2), SEMPO_MUT(1:endday), ':','color', [229, 61, 45]./255, 'linewidth', 2)
ylabel('d''')
xlabel('Day')
ylim([-1 6])
legend(['Wt ',num2str(length(WTtouseold))], ['APP+ ',num2str(length(MUTtouseold))])
xlim([0 19])
title('Old')


%% Memory plots
aw=nanmean(dpblockmat1_Wtyoung,2);
am=nanmean(dpblockmat1_Mutyoung,2);
aw2=nanmean(dpblockmat2_Wtyoung,2);
am2=nanmean(dpblockmat2_Mutyoung,2);
aw3=nanmean(dpblockmat3_Wtyoung,2);
am3=nanmean(dpblockmat3_Mutyoung,2);

bw=nanmean(dpblockmat1_Wtold,2);
bm=nanmean(dpblockmat1_Mutold,2);
bw2=nanmean(dpblockmat2_Wtold,2);
bm2=nanmean(dpblockmat2_Mutold,2);
bw3=nanmean(dpblockmat3_Wtold,2);
bm3=nanmean(dpblockmat3_Mutold,2);

awmt=[aw aw2 aw3];
awmtmm=[am am2 am3];

bwmt=[bw bw2 bw3];
bwmtmm=[bm bm2 bm3];

figure
subplot(2,1,1)
hold all
for wii=1:21
    plot([wii-0.5 wii wii+0.5],awmt(wii,:), 'k-','LineWidth', 2);
    plot([wii-0.5 wii wii+0.5],awmtmm(wii,:), 'r-','LineWidth', 2);
end
xlim([0 22])
% ylim([-1 5])
ylabel('d''')
xlabel('Training day')
title('Young')

subplot(2,1,2)
hold all
for wii=1:21
    plot([wii-0.5 wii wii+0.5],bwmt(wii,:), 'k-','LineWidth', 2);
    plot([wii-0.5 wii wii+0.5],bwmtmm(wii,:), 'r-','LineWidth', 2);
end
xlim([0 22])
% ylim([-1 5])
ylabel('d''')
xlabel('Training day')
title('Old')


%% MEmory plots
colorsml=[108 63 255; 46 255 232; 16 29 66; 246 16 103; 231 98 5; 58 227 140; 255 200 87; 52 101 51; 197 40 61; 35 46 209; 17 153 0; 90 70 146]./255;
dpblockmat1_Wtyoung=NaN(40, length(WTtouseyoung));
dpblockmat2_Wtyoung=NaN(40, length(WTtouseyoung));
dpblockmat3_Wtyoung=NaN(40, length(WTtouseyoung));
dpblockmat1_Mutyoung=NaN(40, length(MUTtouseyoung));
dpblockmat2_Mutyoung=NaN(40, length(MUTtouseyoung));
dpblockmat3_Mutyoung=NaN(40, length(MUTtouseyoung));
dpblockmat1_Wtold=NaN(40, length(WTtouseold));
dpblockmat2_Wtold=NaN(40, length(WTtouseold));
dpblockmat3_Wtold=NaN(40, length(WTtouseold));
dpblockmat1_Mutold=NaN(40, length(MUTtouseold));
dpblockmat2_Mutold=NaN(40, length(MUTtouseold));
dpblockmat3_Mutold=NaN(40, length(MUTtouseold));

figure
subplot(4,1,1)
hold all
i=0;
for ii=WTtouseyoung
    i=i+1;
    dpblockmat1_Wtyoung(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,1);
    dpblockmat2_Wtyoung(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,2);
    dpblockmat3_Wtyoung(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,3);
    for wii=1:length(dp_rblock{ii})
        hold all
        plot([wii-0.4 wii wii+0.4],dp_rblock{ii}(wii,:), '-', 'color', colorsml(i, :),'LineWidth', 2);
    end
end
xlim([0.5 21.5])
ylim([-1 5])
ylabel('d''')
xlabel('Training day')
title('Young Wt')

subplot(4,1,2)
hold all
i=0;
for ii=MUTtouseyoung
    i=i+1;
    dpblockmat1_Mutyoung(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,1);
    dpblockmat2_Mutyoung(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,2);
    dpblockmat3_Mutyoung(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,3);
    for wii=1:length(dp_rblock{ii})
        hold all
        plot([wii-0.4 wii wii+0.4],dp_rblock{ii}(wii,:),  '-', 'color', colorsml(i, :),'LineWidth', 2);
    end
end
xlim([0.5 21.5])
ylim([-1 5])
ylabel('d''')
xlabel('Training day')
title('Young mutants')

subplot(4,1,3)
hold all
i=0;
for ii=WTtouseold
        i=i+1;
    dpblockmat1_Wtold(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,1);
    dpblockmat2_Wtold(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,2);
    dpblockmat3_Wtold(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,3);
    for wii=1:length(dp_rblock{ii})
        hold all
        plot([wii-0.4 wii wii+0.4],dp_rblock{ii}(wii,:),  '-', 'color', colorsml(i, :),'LineWidth', 2);
    end
end
xlim([0.5 21.5])
ylim([-1 5])
ylabel('d''')
xlabel('Training day')
title('Old Wt')

subplot(4,1,4)
hold all
i=0;
for ii=MUTtouseold
    i=i+1;
    dpblockmat1_Mutold(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,1);
    dpblockmat2_Mutold(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,2);
    dpblockmat3_Mutold(1:length(dp_rblock{ii}(:,1)), i)=dp_rblock{ii}(:,3);
    for wii=1:length(dp_rblock{ii})
        hold all
        plot([wii-0.4 wii wii+0.4],dp_rblock{ii}(wii,:),  '-', 'color', colorsml(i, :),'LineWidth', 2);
    end
end
xlim([0.5 21.5])
ylim([-1 5])
ylabel('d''')
xlabel('Training day')
title('Old mutants')

%% Lick training behavior
figure
i=0;
for ii=WTtouseyoung
    i=i+1;
    subplot(4,3,i)
    hold all
    for fi=1:size(LickMatrix{ii},1)
        plot(LickMatrix{ii}(fi,:), 'k-')
        title([num2str(WTtouseyoung(i))])
    end
xlim([1 20])
ylim([0 20])
end

figure
subplot(2,1,1)
hold all
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
for ii=WTtouseyoung
plot(LickMatrix{ii}(5,1:20), 'k-', 'linewidth', 2)
end 
xlim([1 20])
ylim([0 20])
ylabel('Lick #')
xlabel('Bins (1 s)')
subplot(2,1,2)
hold all
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
for ii=MUTtouseyoung
plot(LickMatrix{ii}(5,1:20), 'r-', 'linewidth', 2)
end 
xlim([1 20])
ylim([0 20])
ylabel('Lick #')
xlabel('Bins (1 s)')



i=0;
for ii=WTtouseyoung
    i=i+1;
    for di=[4 5 8 10 15 20]
        if size(LickMatrix{ii},1)>di
            savelickmat{di}(i,:)=LickMatrix{ii}(di,1:20);
        end
    end
end


i=0;
for ii=MUTtouseyoung
    i=i+1;
    for di=[4 5 8 10 15 20]
        if size(LickMatrix{ii},1)>di
            savelickmatmut{di}(i,:)=LickMatrix{ii}(di,1:20);
        end
    end
end

figure
hold all
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
i=0;
for ii=[4 5 8 10 15 20]
    i=i+1;
    subplot(1,6,i)
    hold all
    plot(savelickmatmut{ii}', 'r-', 'linewidth', 2)
    title(['Day ' num2str(ii)])
    xlim([1 20])
    ylim([0 20])
    
end
xlim([1 20])
ylim([0 20])
ylabel('Lick #')
xlabel('Bins (1 s)')

figure
hold all
set(gca,'TickDir','out'); set(gca,'FontSize',12); set(gca,'fontname','arial');
i=0;
for ii=[4 5 8 10 15 20]
    i=i+1;
    subplot(1,6,i)
    hold all
    plot(savelickmat{ii}', 'k-', 'linewidth', 2)
    title(['Day ' num2str(ii)])
    xlim([1 20])
    ylim([0 20])
end

ylabel('Lick #')
xlabel('Bins (1 s)')



