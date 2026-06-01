%% dead period behavior stats 
%first preprocess the data for a given cohort
% cohort=6;
% Behavior_Opto_Dead(cohort)
%% now load the data for each cohort and make one big file called allCohorts
cd('O:\sjk\Behavior\MGBIC_6')
load('deadSummaryData.mat')
% fix the data for 204
% optomeanMat{7,4}(2,1) = optomeanMat{7,16}(5,17);
% optomeanMat{7,5}(2,1) = optomeanMat{7,16}(5,18);
% optomeanMat{7,6}(2,1) = optomeanMat{7,16}(4,15);
% optomeanMat{7,7}(2,1) = optomeanMat{7,16}(4,16);
% optomeanMat(2,:)=[];
delayRates(7,:)=[];
delayRates(12,:)=[];

%%
% save('deadSummaryData.mat','optomeanMat')
% fix the matrix variable
optomeanMat(1,17)={'Matrix as a double'};
for qq=2:size(optomeanMat,1)
    clear newMat
    for rr=1:size(optomeanMat{2,15},1)-4
        if rr==1
            newMat2=optomeanMat{qq,15}{rr+1,1};
            newMat2(:,1)= newMat2(:,1)+rr;
            newMat=vertcat(optomeanMat{qq,15}{rr,1},newMat2);
        else
            newMat2=optomeanMat{qq,15}{rr+1,1};
            newMat2(:,1)= newMat2(:,1)+rr+1;
            newMat=vertcat(newMat,newMat2);
        end
    end
    optomeanMat(qq,17)={newMat};
end
%%
% cohortRange=1:6;
% allCohorts=loadAllOptoCohorts(cohortRange);
SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
% for each animal, record whether it's a test or ctl  % 1 is ctl, 2 is test
condition=[0 2 2 2]; % 198, 203, 204
testIdx=find(condition==2);
allDataTestsOnly=delayRates;

%% make plot to compare percentage correct when light is on vs. off
% Compute percent correct, by session
reinfcolor= [0.4,0.4,0.4];
optocolor=[102/255 178/255 255/255];
firstIdxAnimal=[2,7,12];
for ee=1:length(firstIdxAnimal)
    ii=16;
    allDataTestsOnly{1,ii}='RHit by condition';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,2);
    ii=ii+1;allDataTestsOnly{1,ii}='RFA by condition';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,3);
    ii=ii+1;allDataTestsOnly{1,ii}='Hit Delay 1';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,6);
    ii=ii+1;allDataTestsOnly{1,ii}='FA Delay 1';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,7);
    ii=ii+1;allDataTestsOnly{1,ii}='Hit Delay 2';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,8);
    ii=ii+1;allDataTestsOnly{1,ii}='FA Delay 2';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,9);
    ii=ii+1;allDataTestsOnly{1,ii}='Hit Delay 3';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,10);
    ii=ii+1;allDataTestsOnly{1,ii}='FA Delay 3';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,11);
    ii=ii+1;allDataTestsOnly{1,ii}='Hit Delay 4';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,12);
    ii=ii+1;allDataTestsOnly{1,ii}='FA Delay 4';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,13);
    ii=ii+1;allDataTestsOnly{1,ii}='Hit Delay 5';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,14);
    ii=ii+1;allDataTestsOnly{1,ii}='FA Delay 5';
    allDataTestsOnly{firstIdxAnimal(ee),ii}=allDataTestsOnly(firstIdxAnimal(ee):firstIdxAnimal(ee)+4,15);
end

%%
cd('O:\sjk\Figures\MGB IC Opto');
for jj=1:length(firstIdxAnimal)
    eeFig=figure(jj);hold on;
    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Dead 1
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),18});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),19});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    
    subplot(2,3,1)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    ylabel('percent correct');
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 1']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{firstIdxAnimal(jj),29}=rpc;
    allDataTestsOnly{firstIdxAnimal(jj),30}=opc;

    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 2
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit =cell2mat(allDataTestsOnly{firstIdxAnimal(jj),20});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),21});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,2)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 2']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{firstIdxAnimal(jj),31}=rpc;
    allDataTestsOnly{firstIdxAnimal(jj),32}=opc;

    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 3
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),22});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),23});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    allDataTestsOnly{firstIdxAnimal(jj),33}=rpc;
    allDataTestsOnly{firstIdxAnimal(jj),34}=opc;
    
    subplot(2,3,3)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 3']);
    xticklabels({'light off', 'light on'});

    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 4
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),24});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),25});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    
    subplot(2,3,4)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    ylabel('percent correct');
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 4 ']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{firstIdxAnimal(jj),35}=rpc;
    allDataTestsOnly{firstIdxAnimal(jj),36}=opc;

    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 5
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),26});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),27});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,5)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 5 ']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{firstIdxAnimal(jj),37}=rpc;
    allDataTestsOnly{firstIdxAnimal(jj),38}=opc; 
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{firstIdxAnimal(jj),1}) '_T_MGB_Delay_PercentCorrect_Opto']);
    saveas(gcf,[char(allDataTestsOnly{firstIdxAnimal(jj),1}) '_T_MGB_Delay_PercentCorrect_Opto.png']);
end

close all
%% make plots to compare hit and fa rate with light on vs light off

for jj=1:length(firstIdxAnimal)
    eeFig=figure(jj+3);hold on;
    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Dead 1
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),18});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),19});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,1)
    eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 1']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 2
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit =cell2mat(allDataTestsOnly{firstIdxAnimal(jj),20});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),21});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,2)
    eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 2']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 3
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),22});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),23});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,3)
    eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 3']);
    xticklabels({'hit', 'hit','fa','fa'});
    
    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 4
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),24});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),25});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,4)
    eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' MGB Delay 4']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),16}); % Delay 5
    rfa=cell2mat(allDataTestsOnly{firstIdxAnimal(jj),17});
    ohit = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),26});
    ofa = cell2mat(allDataTestsOnly{firstIdxAnimal(jj),27});
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,5)
    eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{firstIdxAnimal(jj),1} ' Delay 5']);
    xticklabels({'hit', 'hit','fa','fa'});

    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{firstIdxAnimal(jj),1}) '_T_MGB_Delay_HitFARate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{firstIdxAnimal(jj),1}) '_T_MGB_Delay_HitFARate_Opto.png']);
end

%% all animal summary analysis
% group by animal
allDataTestsOnly{1,29}='Delay 1 RPC by Sess';
allDataTestsOnly{1,30}='Delay 1 OPC by Sess';
allDataTestsOnly{1,31}='Delay 2 RPC by Sess';
allDataTestsOnly{1,32}='Delay 2 OPC by Sess';
allDataTestsOnly{1,33}='Delay 3 RPC by Sess';
allDataTestsOnly{1,34}='Delay 3 OPC by Sess';
allDataTestsOnly{1,35}='Delay 4 RPC by Sess';
allDataTestsOnly{1,36}='Delay 4 OPC by Sess';
allDataTestsOnly{1,37}='Delay 5 RPC by Sess';
allDataTestsOnly{1,38}='Delay 5 OPC by Sess';
clear qqq wwFig rpc opc
[allDataTestsOnly]=byAnimalPercentCorrectDelay(allDataTestsOnly,firstIdxAnimal,reinfcolor,optocolor);
close all
clear qqq wwFig rpc opc % now group by session, percent correct for Test
% [allDataTestsOnly]=bySessPercentCorrectDead(allDataTestsOnly,reinfcolor,optocolor);
close all
 %% now do by animal for hit and FA for Test animals
% tempTestsOnly is the IC test animals, mgbTempTestsOnly is the MGB animals
close all
[allDataTestsOnly]=byAnimalHFADead(allDataTestsOnly,firstIdxAnimal,reinfcolor,optocolor);
close all
% bySessHFADead(allDataTestsOnly,firstIdxAnimal,reinfcolor,optocolor);
% close all
%%
% now do the anova
[anovaMat]=anova2OptoPerAnimalDead(allDataTestsOnly,firstIdxAnimal,reinfcolor,optocolor);
anovaInteractionPlotDelay(anovaMat,firstIdxAnimal,allDataTestsOnly) 
% plots the mean and error bars for d'

%% make plot to compare lick latency, for each animal and across animals 
lickLatOptoPerAnimalDead(allDataTestsOnly,reinfcolor,optocolor)

%% AC Expert Data
cd('O:\sjk\Behavior\AC9010Expert'); % this is celine's expert data where she inactivated on 90% of trials
expertACSummaryData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_day.mat'); 
expertACSummaryData=expertACSummaryData.summary_data;
expertACRawData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_lickraster.mat');

plotACExpertData(expertACSummaryData);

%% compare differences in percent correct across each light off vs light on condition
%this is by session
allDataTestsOnly{1,49}='MGB Difference Dead 1';
allDataTestsOnly{1,50}='MGB Difference Dead 2';
allDataTestsOnly{1,51}='MGB Difference Dead 3';
allDataTestsOnly{1,52}='MGB Difference Dead 4';
allDataTestsOnly{1,53}='MGB Difference Dead 5';
% by session
[allDataTestsOnly]=diffPcOptoSessDead(allDataTestsOnly,optocolor);
% by animal
diffPcOptoAnDead(allDataTestsOnly,optocolor);

%% Differences per animal
 % have not updated this yet with the new matrices
ff=2:4;
for ff=2:4
    wwFig=figure(101+ff);
    % subplot(1,3,1)
    qqq=bar([nanmean(allDataTestsOnly{ff,33}) nanmean(allDataTestsOnly{ff,34}) nanmean(allDataTestsOnly{ff,35})]); hold on;
    qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
    scatter(repmat(qqq(1).XEndPoints(1),size(allDataTestsOnly{ff,33},1),1), ...
        allDataTestsOnly{ff,33},'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(qqq(1).XEndPoints(2),size(allDataTestsOnly{ff,34},1),1), ...
        allDataTestsOnly{ff,34},'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(qqq(1).XEndPoints(3),size(allDataTestsOnly{ff,35},1),1), ...
        allDataTestsOnly{ff,35},'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(allDataTestsOnly{ff,33},allDataTestsOnly{ff,34});
    sigstar({[1,2]}, p)
    [h,p,ci,stats] = ttest2(allDataTestsOnly{ff,34},allDataTestsOnly{ff,35});
    sigstar({[2,3]}, p)
    [h,p,ci,stats] = ttest2(allDataTestsOnly{ff,33},allDataTestsOnly{ff,35});
    sigstar({[1,3]}, p)
    title([ char(allDataTestsOnly{ff,1}) 'MGB Inactivation']);
    xticklabels({'Full Trial','Tone','Choice'});
    xlabel('Condition');
    wwFig.Position(3:4)=[325 275];
    ylabel('Light on - Light off');
    saveas(gcf,[ char(allDataTestsOnly{ff,1}) '_Difference_T_MGB_PercentCorrect_Opto']);
    saveas(gcf,[ char(allDataTestsOnly{ff,1}) '_Difference_T_MGB_PercentCorrect_Opto.png']);
end

