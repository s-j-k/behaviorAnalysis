%% behavior stats analysis
%first preprocess the data for a given cohort
cohort=5;
Behavior_Bpod_Opto_SJK(cohort)
%% now load the data for each cohort and make one big file called allCohorts
cohortRange=1:5;
allCohorts=loadAllOptoCohorts(cohortRange);
condition=[0 0 1 0 1 0 1 0 1 ...
    0 2 0 1 0 1 0 2 ... % 164 is a test btu coding as a CTL since little effect/high baseline FA rate
    0 2 0 1 0 1]; % 1 is ctl, 2 is test
testIdx=find(condition==2);
SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;

% for each animal, record whether it's a test or ctl 
condition=[0 0 1 0 1 0 1 0 1 ...
    0 2 0 1 0 1 0 2 ... % 164 is a test btu coding as a CTL since little effect/high baseline FA rate
    0 2 0 1 0 1 0 1 0 1 0 1 0 ... 
    1 0 1 0 1 0 2 0 2 0 2 0]; % cohort 5 last row

testIdx=find(condition==2);
allDataTestsOnly=allCohorts(1,:);
allDataTestsOnly(2:length(testIdx)+1,:)=allCohorts(testIdx,:);
ctlIdx=[21,25,27,31,33,35]; % not all controls are good
allDataCtlOnly=allCohorts(1,:);
allDataCtlOnly(2:length(ctlIdx)+1,:)=allCohorts(ctlIdx,:);

%% plot catch trials
% inherited from celine code... 
% coded as conffa (catch light on false alarm)
% within the rates variable, conffa is the 10th column
% catch trials happen across the last 3 days (1st MGB, then IC, then LOB)
reinfcolor= [0.4,0.4,0.4];
optocolor=[102/255 178/255 255/255];
animalIdx=[3,5,9];
for ww=1:length(animalIdx)
    rates=optomeanMat{animalIdx(ww),27}; % which animal are we plotting
    uu=2;yyyFig=figure('Position', [100 100 800 200]);
    subplot(1,3,1);hold on;
    end1=size(rates,1);
    yyy=bar([rates(end1-uu,1) rates(end1-uu,5); ...
        rates(end1-uu,2) rates(end1-uu,6); 0 rates(end1-uu,10)]);
    % rates(i,:) = [rhit rfa phit pfa ohit ofa ...
    %             rfhit rffa cofffa conffa phit2 pfa2 ...
    %             othit otfa ochit ocfa];
    yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
    yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
    yyy(2).CData=[optocolor;optocolor;optocolor];box on; 
    xticklabels({'hit','','fa','','catch fa'}); legend('light off','light on');
    title([string(optomeanMat{animalIdx(ww),1}) 'MGB Catch']);
    
    uu=1;subplot(1,3,2);
    end1=size(rates,1);
    yyy=bar([rates(end1-uu,1) rates(end1-uu,5); ...
        rates(end1-uu,2) rates(end1-uu,6); 0 rates(end1-uu,10)]);
    yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
    yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
    yyy(2).CData=[optocolor;optocolor;optocolor];
    xticklabels({'hit','fa','catch fa'}); legend('light off','light on');
    title('IC Catch');
        
    subplot(1,3,3);
    end1=size(rates,1);
    yyy=bar([rates(end1,1) rates(end1,5); ...
        rates(end1,2) rates(end1,6); 0 rates(end1,10)]);
    yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
    yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
    yyy(2).CData=[optocolor;optocolor;optocolor];
    xticklabels({'hit','fa','catch fa'}); legend('light off','light on');
    title('LOB Catch');       
    animalName=optomeanMat{animalIdx(ww),1};
        saveas(yyyFig,[char(optomeanMat{animalIdx(ww),1}) '-CatchBar.fig']);
        saveas(yyyFig,[char(optomeanMat{animalIdx(ww),1}) '-CatchBar.png']);
        close(yyyFig);

end

%% make plot to compare percentage correct when light is on vs. off
% Compute percent correct, by session
allDataTestsOnly{1,27}='RPC MGB Full Trial';
allDataTestsOnly{1,28}='OPC MGB Full Trial';
allDataTestsOnly{1,29}='RPC MGB Tone Trial';
allDataTestsOnly{1,30}='OPC MGB Tone Trial';
allDataTestsOnly{1,31}='RPC MGB Choice Trial';
allDataTestsOnly{1,32}='OPC MGB Choice Trial';
allDataTestsOnly{1,33}='RPC IC Full Trial';
allDataTestsOnly{1,34}='OPC IC Full Trial';
allDataTestsOnly{1,35}='RPC IC Tone Trial';
allDataTestsOnly{1,36}='OPC IC Tone Trial';
allDataTestsOnly{1,37}='RPC IC Choice Trial';
allDataTestsOnly{1,38}='OPC IC Choice Trial';
for jj=2:size(allDataTestsOnly,1)
    eeFig=figure(jj);hold on;
    rhit=allDataTestsOnly{jj,10}; % full trial MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,12};
    ofa = allDataTestsOnly{jj,13};
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
    title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'light off', 'light on'});
    allDataTestsOnly{jj,27}=rpc;
    allDataTestsOnly{jj,28}=opc;

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

clear qqq wwFig rpc opc % now group by session, percent correct for Test
[tempTestsOnly,mgbTempTestsOnly,allDataTestsOnly,allDataCtlOnly]=bySessPercentCorrect(allDataTestsOnly,allDataCtlOnly,reinfcolor,optocolor);

 %% now do by animal for hit and FA for Test animals
% tempTestsOnly is the IC test animals, mgbTempTestsOnly is the MGB animals
byAnimalHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,tempTestsOnly,reinfcolor,optocolor);

close all
bySessHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,tempTestsOnly,reinfcolor,optocolor);
%%
% now do the anova....
[aovMGB,statsMGB,aovIC,statsIC]=anova2OptoPerAnimal(mgbTempTestsOnly,tempTestsOnly);
close all
aovMGB
aovIC

%% compare differences in percent correct across each light off vs light on condition
% TODO: needs to be an anova!!!
allDataTestsOnly{1,33}='Difference PC Full Trial';
allDataTestsOnly{1,34}='Difference PC Tone';
allDataTestsOnly{1,35}='Difference PC Choice';
jj=4;
for jj=2:4
    allDataTestsOnly{jj,33}=allDataTestsOnly{jj,28}-allDataTestsOnly{jj,27};
    allDataTestsOnly{jj,34}=allDataTestsOnly{jj,30}-allDataTestsOnly{jj,29};
    allDataTestsOnly{jj,35}=allDataTestsOnly{jj,32}-allDataTestsOnly{jj,31};
end
allDataTestsOnly{4,33}=allDataTestsOnly{4,33}';

diffFull=NaN; diffTone=NaN; diffChoice=NaN;
for jj=2:size(allDataTestsOnly,1)
    diffFull=cat(1,diffFull,allDataTestsOnly{jj,33});
    diffTone=cat(1,diffTone,allDataTestsOnly{jj,34});
    diffChoice=cat(1,diffChoice,allDataTestsOnly{jj,35});
end
wwFig=figure(100);
% subplot(1,3,1)
qqq=bar([nanmean(diffFull) nanmean(diffTone) nanmean(diffChoice)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffFull,1),1), ...
    diffFull,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffTone,1),1), ...
    diffTone,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffChoice,1),1), ...
    diffChoice,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(diffFull,diffTone);
sigstar({[1,2]}, p)
[h,p,ci,stats] = ttest2(diffTone,diffChoice);
sigstar({[2,3]}, p)
[h,p,ci,stats] = ttest2(diffFull,diffChoice);
sigstar({[1,3]}, p)
title(['MGB Inactivation']);
xticklabels({'Full Trial','Tone','Choice'});
xlabel('Condition');
wwFig.Position(3:4)=[325 275];
ylabel('Light on - Light off');
saveas(gcf,['BySess_Difference_T_MGB_PercentCorrect_Opto']);
saveas(gcf,['BySess_Difference_T_MGB_PercentCorrect_Opto.png']);

% by animal

diffFull=NaN; diffTone=NaN; diffChoice=NaN;
for jj=2:size(allDataTestsOnly,1)
    diffFull=cat(1,diffFull,nanmean(allDataTestsOnly{jj,33}));
    diffTone=cat(1,diffTone,nanmean(allDataTestsOnly{jj,34}));
    diffChoice=cat(1,diffChoice,nanmean(allDataTestsOnly{jj,35}));
end
wwFig=figure(101);
% subplot(1,3,1)
qqq=bar([nanmean(diffFull) nanmean(diffTone) nanmean(diffChoice)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[optocolor;optocolor;optocolor;];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(diffFull,1),1), ...
    diffFull,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(2),size(diffTone,1),1), ...
    diffTone,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qqq(1).XEndPoints(3),size(diffChoice,1),1), ...
    diffChoice,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(diffFull,diffTone);
sigstar({[1,2]}, p)
[h,p,ci,stats] = ttest2(diffTone,diffChoice);
sigstar({[2,3]}, p)
[h,p,ci,stats] = ttest2(diffFull,diffChoice);
sigstar({[1,3]}, p)
title(['MGB Inactivation by Animal']);
xticklabels({'Full Trial','Tone','Choice'});
xlabel('Condition');
wwFig.Position(3:4)=[325 275];
ylabel('Light on - Light off');
saveas(gcf,['ByAn_Difference_T_MGB_PercentCorrect_Opto']);
saveas(gcf,['ByAn_Difference_T_MGB_PercentCorrect_Opto.png']);
%%
% per animal

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


%% make plot to compare lick latency 

SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;

% use the raw data 
% tests
nbsubj=2;
for nbsubj=2:size(allDataTestsOnly,1)
    matrix=allDataTestsOnly{nbsubj,26};
    for i=1:max(matrix(:,SESS))
        mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
        sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
        mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        mo_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        so_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        mto_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
        sto_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
        mco_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
        sco_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
        
        mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        mo_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        so_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        mto_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
        sto_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
        mco_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
        sco_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
        
        mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        mlto_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        slto_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        mlco_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        slco_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
                        
        mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        mlto_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        slto_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        mlco_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        slco_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));        
        
    end

    
    switch nbsubj
        case 2
            mgbDays = [1 1 0 0 1 1 0 0 1 1 1 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 1 1 0 0 1 1 0 0 0 1 1];icDays=logical(icDays);
            expRange=1:8; 
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
        case 3
            expRange=22:32;
            mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 0 1 1 0 0 ...
                0 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                1 1 1 0 0 1 1 0 0 1 1 ...
                0 0 0];icDays=logical(icDays);
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
        case 4
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ...
                1 0 1 1 0 1 1 0 1 0 0 1];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ... 
                0 1 0 0 1 0 0 1 0 1 1 0];
            icDays=logical(icDays);
            expRange=20:length(mgbDays); 
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBmr_hit(MGBmr_hit==0)=NaN;
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);           
    end
end

jj=2;     
for jj=2:4
    lickFig=figure(jj+103); 
    subplot(1,3,1); % full trial MGB
    lickbar=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmo_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmo_fa(jj,:))]); hold on;
    error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBso_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBso_fa(jj,:))];
    lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbar)
        xtips=lickbar(kk).XEndPoints;
        ytips=lickbar(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    
    [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmo_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmo_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(1,3,2) % tone MGB
    lickbarT=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmto_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmto_fa(jj,:))]); hold on;
    error=[nanmean(MGBsr_hit(jj,:)) mean(MGBsto_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsto_fa(jj,:))];
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    for kk=1:numel(lickbarT)
        xtips=lickbarT(kk).XEndPoints;
        ytips=lickbarT(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmto_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmto_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(1,3,3); % choice MGB
    lickbarC=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmco_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmco_fa(jj,:))]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBsco_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsco_fa(jj,:))];
    for kk=1:numel(lickbarC)
        xtips=lickbarC(kk).XEndPoints;
        ytips=lickbarC(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmco_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmco_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});
    lickFig.Position(3:4)=[725 275];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickLat_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickLat_Opto.png']);
    
    % LATENCY
    lickFigl=figure(jj+13); subplot(1,3,1); % full trial MGB
    lickbarL=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlo_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlo_fa(jj,:))]); hold on;
    lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslo_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslo_fa(jj,:))];
    for kk=1:numel(lickbarL)
        xtips=lickbarL(kk).XEndPoints;
        ytips=lickbarL(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlo_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlo_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(1,3,2) % tone MGB
    lickbarT=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlto_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlto_fa(jj,:))]); hold on;
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslto_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslto_fa(jj,:))];
    for kk=1:numel(lickbarT)
        xtips=lickbarT(kk).XEndPoints;
        ytips=lickbarT(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlto_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlto_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(1,3,3); % choice MGB
    lickbarC=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlco_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlco_fa(jj,:))]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslco_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslco_fa(jj,:))];
    for kk=1:numel(lickbarC)
        xtips=lickbarC(kk).XEndPoints;
        ytips=lickbarC(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlco_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlco_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});
    lickFigl.Position(3:4)=[725 275];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickRate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_LickRate_Opto.png']);    
end

%% lick rate and latency by animal
lickFig=figure(jj+103); 
for jj=2:size(allDataTestsOnly,1)
    hitLickRate(jj)=nanmean(MGBmr_hit(jj,:));
    ohitLickRate(jj)=nanmean(MGBmo_hit(jj,:));
    faLickRate(jj)=nanmean(MGBmr_fa(jj,:));
    ofaLickRate(jj)=nanmean(MGBmo_fa(jj,:));
    othitLickRate(jj)=nanmean(MGBmto_hit(jj,:));
    otfaLickRate(jj)=nanmean(MGBmto_fa(jj,:));
    ochitLickRate(jj)=nanmean(MGBmco_hit(jj,:));
    ocfaLickRate(jj)=nanmean(MGBmco_fa(jj,:));
end

subplot(1,3,1); % full trial MGB
lickbar=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBso_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBso_fa(jj,:))];
lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbar)
    xtips=lickbar(kk).XEndPoints;
    ytips=lickbar(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(1,3,2) % tone MGB
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
error=[nanmean(MGBsr_hit(jj,:)) mean(MGBsto_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsto_fa(jj,:))];
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(1,3,3); % choice MGB
lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBsco_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsco_fa(jj,:))];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});
lickFig.Position(3:4)=[725 275];
saveas(gcf,['ByAnimal_T_MGB_LickLat_Opto']);
saveas(gcf,['ByAnimal_T_MGB_LickLat_Opto.png']);

% LATENCY
lickFigl=figure(jj+13); subplot(1,3,1); % full trial MGB
lickbarL=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlo_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlo_fa(jj,:))]); hold on;
lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslo_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslo_fa(jj,:))];
for kk=1:numel(lickbarL)
    xtips=lickbarL(kk).XEndPoints;
    ytips=lickbarL(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlo_hit(jj,:));
[h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlo_fa(jj,:));
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(1,3,2) % tone MGB
lickbarT=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlto_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlto_fa(jj,:))]); hold on;
lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslto_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslto_fa(jj,:))];
for kk=1:numel(lickbarT)
    xtips=lickbarT(kk).XEndPoints;
    ytips=lickbarT(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlto_hit(jj,:));
[h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlto_fa(jj,:));
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(1,3,3); % choice MGB
lickbarC=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlco_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlco_fa(jj,:))]); hold on;
lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslco_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslco_fa(jj,:))];
for kk=1:numel(lickbarC)
    xtips=lickbarC(kk).XEndPoints;
    ytips=lickbarC(kk).YEndPoints;
    errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
end
[h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlco_hit(jj,:));
[h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlco_fa(jj,:));
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick rate');
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});
lickFigl.Position(3:4)=[725 275];
saveas(gcf,['ByAnimal_T_MGB_LickRate_Opto']);
saveas(gcf,['ByAnimal_T_MGB_LickRate_Opto.png']);    

