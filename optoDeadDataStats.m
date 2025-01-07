%% dead period behavior stats 
%first preprocess the data for a given cohort
cohort=6;
Behavior_Opto_Dead(cohort)
%% now load the data for each cohort and make one big file called allCohorts
cd('O:\sjk\Behavior\MGBIC_6')
load('deadSummaryData.mat')
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
    allDataTestsOnly{jj,17}=rpc;
    allDataTestsOnly{jj,18}=opc;

    rhit=allDataTestsOnly{jj,10}; % tone MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,14};
    ofa = allDataTestsOnly{jj,15};
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
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,29}=rpc;
    allDataTestsOnly{jj,30}=opc;

    rhit=allDataTestsOnly{jj,10}; % choice MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,16};
    ofa = allDataTestsOnly{jj,17};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    allDataTestsOnly{jj,31}=rpc;
    allDataTestsOnly{jj,32}=opc;
    
    subplot(2,3,3)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataTestsOnly{jj,1} ' MGB Choice Inactivation']);
    xticklabels({'light off', 'light on'});

    % now do the same for the IC
    rhit=allDataTestsOnly{jj,18}; % full trial IC
    rfa=allDataTestsOnly{jj,19};
    ohit = allDataTestsOnly{jj,20};
    ofa = allDataTestsOnly{jj,21};
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
    title([allDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,33}=rpc;
    allDataTestsOnly{jj,34}=opc;

    rhit=allDataTestsOnly{jj,18}; % tone IC
    rfa=allDataTestsOnly{jj,19};
    ohit = allDataTestsOnly{jj,22};
    ofa = allDataTestsOnly{jj,23};
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
    title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,35}=rpc;
    allDataTestsOnly{jj,36}=opc;

    rhit=allDataTestsOnly{jj,18}; % choice IC
    rfa=allDataTestsOnly{jj,19};
    ohit = allDataTestsOnly{jj,24};
    ofa = allDataTestsOnly{jj,25};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    allDataTestsOnly{jj,37}=rpc;
    allDataTestsOnly{jj,38}=opc;
    
    subplot(2,3,6)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataTestsOnly{jj,1} ' IC Choice Inactivation']);
    xticklabels({'light off', 'light on'});    
    
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_PercentCorrect_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_PercentCorrect_Opto.png']);
    
end

allDataCtlOnly{1,27}='RPC MGB Full Trial';
allDataCtlOnly{1,28}='OPC MGB Full Trial';
allDataCtlOnly{1,29}='RPC MGB Tone Trial';
allDataCtlOnly{1,30}='OPC MGB Tone Trial';
allDataCtlOnly{1,31}='RPC MGB Choice Trial';
allDataCtlOnly{1,32}='OPC MGB Choice Trial';
allDataCtlOnly{1,33}='RPC IC Full Trial';
allDataCtlOnly{1,34}='OPC IC Full Trial';
allDataCtlOnly{1,35}='RPC IC Tone Trial';
allDataCtlOnly{1,36}='OPC IC Tone Trial';
allDataCtlOnly{1,37}='RPC IC Choice Trial';
allDataCtlOnly{1,38}='OPC IC Choice Trial';
for jj=2:size(allDataCtlOnly,1) % control animals
    eeFig=figure(jj+10);hold on;
    rhit=allDataCtlOnly{jj,10}; % full trial MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,12};
    ofa = allDataCtlOnly{jj,13};
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
    title([allDataCtlOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataCtlOnly{jj,27}=rpc;
    allDataCtlOnly{jj,28}=opc;

    rhit=allDataCtlOnly{jj,10}; % tone MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,14};
    ofa = allDataCtlOnly{jj,15};
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
    title([allDataCtlOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataCtlOnly{jj,29}=rpc;
    allDataCtlOnly{jj,30}=opc;

    rhit=allDataCtlOnly{jj,10}; % choice MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,16};
    ofa = allDataCtlOnly{jj,17};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    allDataCtlOnly{jj,31}=rpc;
    allDataCtlOnly{jj,32}=opc;
    subplot(2,3,3)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataCtlOnly{jj,1} ' MGB Choice Inactivation']);
    xticklabels({'light off', 'light on'});
    
    rhit=allDataCtlOnly{jj,18}; % full IC
    rfa=allDataCtlOnly{jj,19};
    ohit = allDataCtlOnly{jj,20};
    ofa = allDataCtlOnly{jj,21};
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
    title([allDataCtlOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataCtlOnly{jj,33}=rpc;
    allDataCtlOnly{jj,34}=opc;

    rhit=allDataCtlOnly{jj,18}; % tone IC
    rfa=allDataCtlOnly{jj,19};
    ohit = allDataCtlOnly{jj,22};
    ofa = allDataCtlOnly{jj,23};
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
    title([allDataCtlOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataCtlOnly{jj,35}=rpc;
    allDataCtlOnly{jj,36}=opc;

    rhit=allDataCtlOnly{jj,18}; % choice IC
    rfa=allDataCtlOnly{jj,19};
    ohit = allDataCtlOnly{jj,24};
    ofa = allDataCtlOnly{jj,25};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    allDataCtlOnly{jj,37}=rpc;
    allDataCtlOnly{jj,38}=opc;
    
    subplot(2,3,6)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([allDataCtlOnly{jj,1} ' IC Choice Inactivation']);
    xticklabels({'light off', 'light on'});
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_IC_PercentCorrect_Opto']);
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_IC_PercentCorrect_Opto.png']);
end
close all
%% make plots to compare hit and fa rate with light on vs light off

for jj=2:size(allDataTestsOnly,1)
    eeFig=figure(jj+3);hold on;
    rhit=allDataTestsOnly{jj,10}; % full trial MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,12};
    ofa = allDataTestsOnly{jj,13};
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
    title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataTestsOnly{jj,10}; % tone MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,14};
    ofa = allDataTestsOnly{jj,15};
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
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataTestsOnly{jj,10}; % choice MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,16};
    ofa = allDataTestsOnly{jj,17};
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
    title([allDataTestsOnly{jj,1} ' MGB Choice Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});
    
    rhit=allDataTestsOnly{jj,18}; % full trial IC
    rfa=allDataTestsOnly{jj,19};
    ohit = allDataTestsOnly{jj,20};
    ofa = allDataTestsOnly{jj,21};
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
    title([allDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataTestsOnly{jj,18}; % tone IC
    rfa=allDataTestsOnly{jj,19};
    ohit = allDataTestsOnly{jj,22};
    ofa = allDataTestsOnly{jj,23};
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
    title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataTestsOnly{jj,18}; % choice IC
    rfa=allDataTestsOnly{jj,19};
    ohit = allDataTestsOnly{jj,24};
    ofa = allDataTestsOnly{jj,25};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,6)
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
    title([allDataTestsOnly{jj,1} ' IC Choice Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_HitFARate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_HitFARate_Opto.png']);
end

clear eee eeFig 
%CONTROLS
for jj=2:size(allDataCtlOnly,1)
    eeFig=figure(jj+10);hold on;
    rhit=allDataCtlOnly{jj,10}; % full trial MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,12};
    ofa = allDataCtlOnly{jj,13};
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
    title([allDataCtlOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataCtlOnly{jj,10}; % tone MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,14};
    ofa = allDataCtlOnly{jj,15};
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
    title([allDataCtlOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataCtlOnly{jj,10}; % choice MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,16};
    ofa = allDataCtlOnly{jj,17};
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
    title([allDataCtlOnly{jj,1} ' MGB Choice Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});
    
    rhit=allDataCtlOnly{jj,18}; % full trial IC
    rfa=allDataCtlOnly{jj,19};
    ohit = allDataCtlOnly{jj,20};
    ofa = allDataCtlOnly{jj,21};
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
    title([allDataCtlOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataCtlOnly{jj,18}; % toneIC
    rfa=allDataCtlOnly{jj,19};
    ohit = allDataCtlOnly{jj,22};
    ofa = allDataCtlOnly{jj,23};
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
    title([allDataCtlOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    rhit=allDataCtlOnly{jj,18}; % choice IC
    rfa=allDataCtlOnly{jj,19};
    ohit = allDataCtlOnly{jj,24};
    ofa = allDataCtlOnly{jj,25};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,6)
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
    title([allDataCtlOnly{jj,1} ' IC Choice Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});   
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_IC_HitFARate_Opto']);
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_IC_HitFARate_Opto.png']);
end

%% all animal summary analysis
% group by animal

clear qqq wwFig rpc opc
[tempTestsOnly,mgbTempTestsOnly,allDataTestsOnly,allDataCtlOnly]=byAnimalPercentCorrect(allDataTestsOnly,allDataCtlOnly,reinfcolor,optocolor);
close all
clear qqq wwFig rpc opc % now group by session, percent correct for Test
[allDataTestsOnly,allDataCtlOnly]=bySessPercentCorrect(allDataTestsOnly,allDataCtlOnly,tempTestsOnly,mgbTempTestsOnly,reinfcolor,optocolor);
close all
 %% now do by animal for hit and FA for Test animals
% tempTestsOnly is the IC test animals, mgbTempTestsOnly is the MGB animals
byAnimalHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,tempTestsOnly,reinfcolor,optocolor);
close all

bySessHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,tempTestsOnly,reinfcolor,optocolor);
close all
%%
% now do the anova....
[aovMGB,statsMGB,aovIC,statsIC]=anova2OptoPerAnimal(mgbTempTestsOnly,tempTestsOnly);
% close all
aovMGB
% aovIC

%% compare differences in percent correct across each light off vs light on condition
%this is by session
allDataTestsOnly{1,39}='MGB Difference PC Full Trial';
allDataTestsOnly{1,40}='MGB Difference PC Tone';
allDataTestsOnly{1,41}='MGB Difference PC Choice';
allDataTestsOnly{1,42}='IC Difference PC Full Trial';
allDataTestsOnly{1,43}='IC Difference PC Tone';
allDataTestsOnly{1,44}='IC Difference PC Choice';
% by session
[allDataTestsOnly,mgbTempTestsOnly,tempTestsOnly]=diffPcOptoSess(allDataTestsOnly,mgbTempTestsOnly,tempTestsOnly,optocolor);
% by animal
diffPcOptoAn(allDataTestsOnly,optocolor);

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
lickLatOptoPerAnimal(mgbTempTestsOnly,tempTestsOnly,allDataTestsOnly,reinfcolor,optocolor)

%% AC Expert Data
cd('O:\sjk\Behavior\AC9010Expert'); % this is celine's expert data where she inactivated on 90% of trials
expertACSummaryData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_day.mat'); 
expertACSummaryData=expertACSummaryData.summary_data;
expertACRawData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_lickraster.mat');

plotACExpertData(expertACSummaryData);