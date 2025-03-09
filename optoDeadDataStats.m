%% dead period behavior stats 
%first preprocess the data for a given cohort
cohort=6;
Behavior_Opto_Dead(cohort)
%% now load the data for each cohort and make one big file called allCohorts
cd('O:\sjk\Behavior\MGBIC_6')
load('deadSummaryData.mat')
% fix the data for 204
% optomeanMat{7,4}(2,1) = optomeanMat{7,16}(5,17);
% optomeanMat{7,5}(2,1) = optomeanMat{7,16}(5,18);
% optomeanMat{7,6}(2,1) = optomeanMat{7,16}(4,15);
% optomeanMat{7,7}(2,1) = optomeanMat{7,16}(4,16);
% optomeanMat(2,:)=[];
% optomeanMat(3,:)=[];
% optomeanMat(4,:)=[];
% save('deadSummaryData.mat','optomeanMat')
%%
% cohortRange=1:6;
% allCohorts=loadAllOptoCohorts(cohortRange);
SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
% for each animal, record whether it's a test or ctl  % 1 is ctl, 2 is test
condition=[0 2 2 2]; % 198, 203, 204
testIdx=find(condition==2);
allDataTestsOnly=optomeanMat;

%% make plot to compare percentage correct when light is on vs. off
% Compute percent correct, by session
reinfcolor= [0.4,0.4,0.4];
optocolor=[102/255 178/255 255/255];
allDataTestsOnly{1,17}='RHit by Condition';
allDataTestsOnly{1,18}='RFA by Condition';
allDataTestsOnly{1,19}='RPC MGB Dead 1';
allDataTestsOnly{1,20}='OPC MGB Dead 1';
allDataTestsOnly{1,21}='RPC MGB Dead 2';
allDataTestsOnly{1,22}='OPC MGB Dead 2';
allDataTestsOnly{1,23}='RPC MGB Dead 3';
allDataTestsOnly{1,24}='OPC MGB Dead 3';
allDataTestsOnly{1,25}='RPC MGB Dead 4';
allDataTestsOnly{1,26}='OPC MGB Dead 4';
allDataTestsOnly{1,27}='RPC MGB Dead 5';
allDataTestsOnly{1,28}='OPC MGB Dead 5';
% organize the reinforced hit sessions by the dead periods (dead 1 through
% dead 5 is each row)
for jj = 2:size(allDataTestsOnly,1)
    allDataTestsOnly{jj,17}(1,1)=allDataTestsOnly{jj,2}(1,1);
    allDataTestsOnly{jj,17}(1,2)=allDataTestsOnly{jj,2}(4,1);
    allDataTestsOnly{jj,17}(2,1)=allDataTestsOnly{jj,2}(1,1);
    allDataTestsOnly{jj,17}(2,2)=allDataTestsOnly{jj,2}(5,1);
    allDataTestsOnly{jj,17}(3,1)=allDataTestsOnly{jj,2}(2,1);
    allDataTestsOnly{jj,17}(3,2)=allDataTestsOnly{jj,2}(3,1);
    allDataTestsOnly{jj,17}(4,1)=allDataTestsOnly{jj,2}(2,1);
    allDataTestsOnly{jj,17}(4,2)=allDataTestsOnly{jj,2}(3,1);
    allDataTestsOnly{jj,17}(5,1)=allDataTestsOnly{jj,2}(4,1);
    allDataTestsOnly{jj,17}(5,2)=allDataTestsOnly{jj,2}(5,1);
    
    allDataTestsOnly{jj,18}(1,1)=allDataTestsOnly{jj,3}(1,1);
    allDataTestsOnly{jj,18}(1,2)=allDataTestsOnly{jj,3}(4,1);
    allDataTestsOnly{jj,18}(2,1)=allDataTestsOnly{jj,3}(1,1);
    allDataTestsOnly{jj,18}(2,2)=allDataTestsOnly{jj,3}(5,1);
    allDataTestsOnly{jj,18}(3,1)=allDataTestsOnly{jj,3}(2,1);
    allDataTestsOnly{jj,18}(3,2)=allDataTestsOnly{jj,3}(3,1);
    allDataTestsOnly{jj,18}(4,1)=allDataTestsOnly{jj,3}(2,1);
    allDataTestsOnly{jj,18}(4,2)=allDataTestsOnly{jj,3}(3,1);
    allDataTestsOnly{jj,18}(5,1)=allDataTestsOnly{jj,3}(4,1);
    allDataTestsOnly{jj,18}(5,2)=allDataTestsOnly{jj,3}(5,1);
end

%%
cd('O:\sjk\Figures\MGB IC Opto');
for jj=2:size(allDataTestsOnly,1)
    eeFig=figure(jj);hold on;
    rhit=allDataTestsOnly{jj,17}(1,:); % Dead 1
    rfa=allDataTestsOnly{jj,18}(1,:);
    ohit = allDataTestsOnly{jj,4};
    ofa = allDataTestsOnly{jj,5};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 1 Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,19}=rpc;
    allDataTestsOnly{jj,20}=opc;

    rhit=allDataTestsOnly{jj,17}(2,:); % Dead 2
    rfa=allDataTestsOnly{jj,18}(2,:);
    ohit = allDataTestsOnly{jj,6};
    ofa = allDataTestsOnly{jj,7};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 2 Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,21}=rpc;
    allDataTestsOnly{jj,22}=opc;

    rhit=allDataTestsOnly{jj,17}(3,:); % Dead 3
    rfa=allDataTestsOnly{jj,18}(3,:);
    ohit = allDataTestsOnly{jj,8};
    ofa = allDataTestsOnly{jj,9};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    allDataTestsOnly{jj,23}=rpc;
    allDataTestsOnly{jj,24}=opc;
    
    subplot(2,3,3)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataTestsOnly{jj,1} ' MGB Dead 3 Inactivation']);
    xticklabels({'light off', 'light on'});

    rhit=allDataTestsOnly{jj,17}(4,:); % Dead 4
    rfa=allDataTestsOnly{jj,18}(4,:);
    ohit = allDataTestsOnly{jj,10};
    ofa = allDataTestsOnly{jj,11};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 4 Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,25}=rpc;
    allDataTestsOnly{jj,26}=opc;

    rhit=allDataTestsOnly{jj,17}(5,:); % Dead 5
    rfa=allDataTestsOnly{jj,18}(5,:);
    ohit = allDataTestsOnly{jj,12};
    ofa = allDataTestsOnly{jj,13};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 5 Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,27}=rpc;
    allDataTestsOnly{jj,28}=opc; 
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_Dead_PercentCorrect_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_Dead_PercentCorrect_Opto.png']);
end

close all
%% make plots to compare hit and fa rate with light on vs light off

for jj=2:size(allDataTestsOnly,1)
    eeFig=figure(jj+3);hold on;
    rhit=allDataTestsOnly{jj,17}(1,:); % Dead 1
    rfa=allDataTestsOnly{jj,18}(1,:);
    ohit = allDataTestsOnly{jj,4};
    ofa = allDataTestsOnly{jj,5};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 1 Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataTestsOnly{jj,17}(2,:); % Dead 2
    rfa=allDataTestsOnly{jj,18}(2,:);
    ohit = allDataTestsOnly{jj,6};
    ofa = allDataTestsOnly{jj,7};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 2 Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataTestsOnly{jj,17}(3,:); % Dead 3
    rfa=allDataTestsOnly{jj,18}(3,:);
    ohit = allDataTestsOnly{jj,8};
    ofa = allDataTestsOnly{jj,9};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 3 Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});
    
    rhit=allDataTestsOnly{jj,17}(4,:); % Dead 4
    rfa=allDataTestsOnly{jj,18}(4,:);
    ohit = allDataTestsOnly{jj,10};
    ofa = allDataTestsOnly{jj,11};
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
    title([allDataTestsOnly{jj,1} ' MGB Dead 4 Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataTestsOnly{jj,17}(5,:); % Dead 5
    rfa=allDataTestsOnly{jj,18}(5,:);
    ohit = allDataTestsOnly{jj,12};
    ofa = allDataTestsOnly{jj,13};
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
    title([allDataTestsOnly{jj,1} ' Dead 5 Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_Dead_HitFARate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_Dead_HitFARate_Opto.png']);
end

%% all animal summary analysis
% group by animal
allDataTestsOnly{1,29}='Dead 1 RPC by Animal';
allDataTestsOnly{1,30}='Dead 1 OPC by Animal';
allDataTestsOnly{1,31}='Dead 2 RPC by Animal';
allDataTestsOnly{1,32}='Dead 2 OPC by Animal';
allDataTestsOnly{1,33}='Dead 3 RPC by Animal';
allDataTestsOnly{1,34}='Dead 3 OPC by Animal';
allDataTestsOnly{1,35}='Dead 4 RPC by Animal';
allDataTestsOnly{1,36}='Dead 4 OPC by Animal';
allDataTestsOnly{1,37}='Dead 5 RPC by Animal';
allDataTestsOnly{1,38}='Dead 5 OPC by Animal';
allDataTestsOnly{1,39}='Dead 1 RPC by Sess';
allDataTestsOnly{1,40}='Dead 1 OPC by Sess';
allDataTestsOnly{1,41}='Dead 2 RPC by Sess';
allDataTestsOnly{1,42}='Dead 2 OPC by Sess';
allDataTestsOnly{1,43}='Dead 3 RPC by Sess';
allDataTestsOnly{1,44}='Dead 3 OPC by Sess';
allDataTestsOnly{1,45}='Dead 4 RPC by Sess';
allDataTestsOnly{1,46}='Dead 4 OPC by Sess';
allDataTestsOnly{1,47}='Dead 5 RPC by Sess';
allDataTestsOnly{1,48}='Dead 5 OPC by Sess';
clear qqq wwFig rpc opc
[allDataTestsOnly]=byAnimalPercentCorrectDead(allDataTestsOnly,reinfcolor,optocolor);
close all
clear qqq wwFig rpc opc % now group by session, percent correct for Test
% [allDataTestsOnly]=bySessPercentCorrectDead(allDataTestsOnly,reinfcolor,optocolor);
close all
 %% now do by animal for hit and FA for Test animals
% tempTestsOnly is the IC test animals, mgbTempTestsOnly is the MGB animals
byAnimalHFADead(allDataTestsOnly,reinfcolor,optocolor);
close all

% bySessHFADead(allDataTestsOnly,reinfcolor,optocolor);
% close all
%%
% now do the anova
[aovMGB,statsMGB]=anova2OptoPerAnimalDead(allDataTestsOnly);
aovMGB

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

%% make plot to compare lick latency, for each animal and across animals 
lickLatOptoPerAnimalDead(allDataTestsOnly,reinfcolor,optocolor)

%% AC Expert Data
cd('O:\sjk\Behavior\AC9010Expert'); % this is celine's expert data where she inactivated on 90% of trials
expertACSummaryData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_day.mat'); 
expertACSummaryData=expertACSummaryData.summary_data;
expertACRawData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_lickraster.mat');

plotACExpertData(expertACSummaryData);