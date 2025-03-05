%% behavior stats analysis
%first preprocess the data for a given cohort
% cohort=5;
% Behavior_Bpod_Opto_SJK(cohort)
%% now load the data for each cohort and make one big file called allCohorts
cohortRange=1:6;
allCohorts=loadAllOptoCohorts(cohortRange);
SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;

% for each animal, record whether it's a test or ctl  % 1 is ctl, 2 is test
condition=[0 0 1 0 1 0 1 0 1 ... %158-162
    0 2 0 1 0 1 0 2 ... % 164 is a test btu coding as a CTL since little effect/high baseline FA rate
    0 2 0 1 0 1 0 1 ...% 176 - 179
    0 1 0 1 ... %183 - 186
    0 1 0 1 0 1 ...
    0 1 0 2 0 2 0 2]; % cohort 6, 198, 203, 204 
%     0 0 0 0 0 0 ...  % 200, 201, 202 who did not learn
% condition=[0 0 0 0 0 0 0 0 0 0 0 0 2]; % cohort 6, 198, 203, 204 

testIdx=find(condition==2);
allDataTestsOnly=allCohorts(1,:);
allDataTestsOnly(2:length(testIdx)+1,:)=allCohorts(testIdx,:);
ctlIdx=[21,25,27,31,33,35]; % not all controls are good
allDataCtlOnly=allCohorts(1,:);
allDataCtlOnly(2:length(ctlIdx)+1,:)=allCohorts(ctlIdx,:);
icDataTestsOnly=allCohorts(1,:);
icDataTestsOnly(2,:)=allCohorts(37,:);
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
allDataTestsOnly{1,39}='stats MGB Full, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,40}='stats MGB tone, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,41}='stats MGB choice, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,42}='stats IC full, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,43}='stats IC tone, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,44}='stats IC choice, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,45}='stats MGB Full dprime, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,46}='stats MGB tone dprime, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,47}='stats MGB choice dprime, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,48}='stats IC Full dprime, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,49}='stats IC Tone dprime, p hit vs opto hit and fa vs opto fa';
allDataTestsOnly{1,50}='stats IC choice dprime, p hit vs opto hit and fa vs opto fa';

icDataTestsOnly{1,27}='RPC MGB Full Trial';
icDataTestsOnly{1,28}='OPC MGB Full Trial';
icDataTestsOnly{1,29}='RPC MGB Tone Trial';
icDataTestsOnly{1,30}='OPC MGB Tone Trial';
icDataTestsOnly{1,31}='RPC MGB Choice Trial';
icDataTestsOnly{1,32}='OPC MGB Choice Trial';
icDataTestsOnly{1,33}='RPC IC Full Trial';
icDataTestsOnly{1,34}='OPC IC Full Trial';
icDataTestsOnly{1,35}='RPC IC Tone Trial';
icDataTestsOnly{1,36}='OPC IC Tone Trial';
icDataTestsOnly{1,37}='RPC IC Choice Trial';
icDataTestsOnly{1,38}='OPC IC Choice Trial';

%% plot catch trials
% inherited from celine code... 
% coded as conffa (catch light on false alarm)
% within the rates variable, conffa is the 10th column
% catch trials happen across the last 3 days (1st MGB, then IC, then LOB)

% reinfcolor= [0.4,0.4,0.4];
% optocolor=[102/255 178/255 255/255];
% animalIdx=[3,5,9];
% for ww=1:length(animalIdx)
%     rates=optomeanMat{animalIdx(ww),27}; % which animal are we plotting
%     uu=2;yyyFig=figure('Position', [100 100 800 200]);
%     subplot(1,3,1);hold on;
%     end1=size(rates,1);
%     yyy=bar([rates(end1-uu,1) rates(end1-uu,5); ...
%         rates(end1-uu,2) rates(end1-uu,6); 0 rates(end1-uu,10)]);
%     % rates(i,:) = [rhit rfa phit pfa ohit ofa ...
%     %             rfhit rffa cofffa conffa phit2 pfa2 ...
%     %             othit otfa ochit ocfa];
%     yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
%     yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
%     yyy(2).CData=[optocolor;optocolor;optocolor];box on; 
%     xticklabels({'hit','','fa','','catch fa'}); legend('light off','light on');
%     title([string(optomeanMat{animalIdx(ww),1}) 'MGB Catch']);
%     
%     uu=1;subplot(1,3,2);
%     end1=size(rates,1);
%     yyy=bar([rates(end1-uu,1) rates(end1-uu,5); ...
%         rates(end1-uu,2) rates(end1-uu,6); 0 rates(end1-uu,10)]);
%     yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
%     yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
%     yyy(2).CData=[optocolor;optocolor;optocolor];
%     xticklabels({'hit','fa','catch fa'}); legend('light off','light on');
%     title('IC Catch');
%         
%     subplot(1,3,3);
%     end1=size(rates,1);
%     yyy=bar([rates(end1,1) rates(end1,5); ...
%         rates(end1,2) rates(end1,6); 0 rates(end1,10)]);
%     yyy(1).FaceColor='flat';yyy(2).FaceColor='flat';
%     yyy(1).CData=[reinfcolor;reinfcolor;reinfcolor];ylim([0 1]);
%     yyy(2).CData=[optocolor;optocolor;optocolor];
%     xticklabels({'hit','fa','catch fa'}); legend('light off','light on');
%     title('LOB Catch');       
%     animalName=optomeanMat{animalIdx(ww),1};
%         saveas(yyyFig,[char(optomeanMat{animalIdx(ww),1}) '-CatchBar.fig']);
%         saveas(yyyFig,[char(optomeanMat{animalIdx(ww),1}) '-CatchBar.png']);
%         close(yyyFig);
% 
% end

%% make plot to compare percentage correct when light is on vs. off
% Compute percent correct, by session
reinfcolor= [0.4,0.4,0.4];
optocolor=[102/255 178/255 255/255];

close all
xVector = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
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
    line([1 2],[rpc opc], 'LineWidth', 0.5, 'Color', [0 0 0]);
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
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_PercentCorrect_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_PercentCorrect_Opto.png']);    
end

for jj=2:size(icDataTestsOnly,1)
    eeFig=figure(jj+20);hold on;
    % now do the same for the IC
    rhit=icDataTestsOnly{jj,18}; % full trial IC
    rfa=icDataTestsOnly{jj,19};
    ohit = icDataTestsOnly{jj,20};
    ofa = icDataTestsOnly{jj,21};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    
    subplot(2,3,4)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6) rpc(7) opc(7) ...
    rpc(8) opc(8)];
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xVector=xVector(~isnan(percentCorrect));
    percentCorrect=percentCorrect(~isnan(percentCorrect));
    line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    ylabel('percent correct');
    title([icDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'light off', 'light on'});
    icDataTestsOnly{jj,33}=rpc;
    icDataTestsOnly{jj,34}=opc;

    rhit=icDataTestsOnly{jj,18}; % tone IC
    rfa=icDataTestsOnly{jj,19};
    ohit = icDataTestsOnly{jj,22};
    ofa = icDataTestsOnly{jj,23};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,5)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xVector=xVector(~isnan(percentCorrect));
    percentCorrect=percentCorrect(~isnan(percentCorrect));
    line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([icDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'light off', 'light on'});
    icDataTestsOnly{jj,35}=rpc;
    icDataTestsOnly{jj,36}=opc;

    rhit=icDataTestsOnly{jj,18}; % choice IC
    rfa=icDataTestsOnly{jj,19};
    ohit = icDataTestsOnly{jj,24};
    ofa = icDataTestsOnly{jj,25};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    icDataTestsOnly{jj,37}=rpc;
    icDataTestsOnly{jj,38}=opc;
    
    subplot(2,3,6)
    eee=bar([nanmean(rpc) nanmean(opc)]); hold on;
    percentCorrect = [rpc(1) opc(1) rpc(2) opc(2) rpc(3) opc(3) ...
    rpc(4) opc(4) rpc(5) opc(5) rpc(6) opc(6)];
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor];hold on;
    scatter(repmat(eee(1).XEndPoints(1),size(rpc,1),1), ...
        rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(opc,1),1), ...
        opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xVector=xVector(~isnan(percentCorrect));
    percentCorrect=percentCorrect(~isnan(percentCorrect));
    line([1 2],percentCorrect(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([1 2],percentCorrect(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,p,ci,stats] = ttest2(rpc,opc);
    sigstar({[1,2]}, p)
    title([icDataTestsOnly{jj,1} ' IC Choice Inactivation']);
    xticklabels({'light off', 'light on'});    
    
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(icDataTestsOnly{jj,1}) '_T_IC_PercentCorrect_Opto']);
    saveas(gcf,[char(icDataTestsOnly{jj,1}) '_T_IC_PercentCorrect_Opto.png']);
    
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
close all
for jj=2:size(allDataTestsOnly,1)
    eeFig=figure(jj+3);hold on;
    rhit=allDataTestsOnly{jj,10}; % full trial MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,12};
    ofa = allDataTestsOnly{jj,13};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,1)
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
    xticklabels({'hit','fa','hit','fa'});

    rhit=allDataTestsOnly{jj,10}; % tone MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,14};
    ofa = allDataTestsOnly{jj,15};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,2)
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
    xticklabels({'hit','fa','hit','fa'});

    rhit=allDataTestsOnly{jj,10}; % choice MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,16};
    ofa = allDataTestsOnly{jj,17};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,3)
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataTestsOnly{jj,1} ' MGB Choice Inactivation']);
    xticklabels({'hit','fa','hit','fa'});
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_HitFARate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_HitFARate_Opto.png']);
    
    
    figure(jj+8);
    rhit=allDataTestsOnly{jj,10}; % full trial MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,12};
    ofa = allDataTestsOnly{jj,13};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    [rdp,c]=dprime_simple(rhit,rfa);
    [odp,c]=dprime_simple(ohit,ofa);
    subplot(1,3,1); plot(rdp,'-o','Color',reinfcolor); hold on;plot(odp,'-o','Color',optocolor);
    title('Full');
    rhit=allDataTestsOnly{jj,10}; % tone MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,14};
    ofa = allDataTestsOnly{jj,15};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    [rdp,c]=dprime_simple(rhit,rfa);
    [odp,c]=dprime_simple(ohit,ofa);
    subplot(1,3,1); plot(rdp,'-o','Color',reinfcolor); hold on;plot(odp,'-o','Color',optocolor);
    title('Stimulus');
    rhit=allDataTestsOnly{jj,10}; % choice MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,16};
    ofa = allDataTestsOnly{jj,17};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    [rdp,c]=dprime_simple(rhit,rfa);
    [odp,c]=dprime_simple(ohit,ofa);
    subplot(1,3,1); plot(rdp,'-o','Color',reinfcolor); hold on;plot(odp,'-o','Color',optocolor);
    title('Choice');
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_dprime_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_dprime_Opto.png']);
end

for jj=2:size(icDataTestsOnly,1)
    eeFig=figure(jj+21);hold on;
    rhit=icDataTestsOnly{jj,18}; % full trial IC
    rfa=icDataTestsOnly{jj,19};
    ohit = icDataTestsOnly{jj,20};
    ofa = icDataTestsOnly{jj,21};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,4)
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('rate');
    title([icDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
    xticklabels({'hit','fa','hit','fa'});

    rhit=icDataTestsOnly{jj,18}; % tone IC
    rfa=icDataTestsOnly{jj,19};
    ohit = icDataTestsOnly{jj,22};
    ofa = icDataTestsOnly{jj,23};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(2,3,5)
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('rate');
    title([icDataTestsOnly{jj,1} ' IC Tone Inactivation']);
    xticklabels({'hit','fa','hit','fa'});

    rhit=icDataTestsOnly{jj,18}; % choice IC
    rfa=icDataTestsOnly{jj,19};
    ohit = icDataTestsOnly{jj,24};
    ofa = icDataTestsOnly{jj,25};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(2,3,6)
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('rate');
    title([icDataTestsOnly{jj,1} ' IC Choice Inactivation']);
    xticklabels({'hit','fa','hit','fa'});
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(icDataTestsOnly{jj,1}) '_T_IC_HitFARate_Opto']);
    saveas(gcf,[char(icDataTestsOnly{jj,1}) '_T_IC_HitFARate_Opto.png']);
end

clear eee eeFig 
close all
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
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
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
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
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
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
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
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
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
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
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
    eee=bar([nanmean(rhit) nanmean(rfa) nanmean(ohit) nanmean(ofa)]); hold on;
    eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
        rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(2),size(rfa,1),2), ...
        rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(eee(1).XEndPoints(3),size(ohit,1),1), ...
        ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
        ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    line([1 2],[rhit rfa], 'LineWidth', 0.5, 'Color', [0 0 0]);
    line([3 4],[ohit ofa],'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(rhit,ohit);
    [h,pFA,ci,stats] = ttest2(rfa,ofa);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('rate');
    title([allDataCtlOnly{jj,1} ' IC Choice Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});   
    
    eeFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_IC_HitFARate_Opto']);
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_IC_HitFARate_Opto.png']);
end
close all
%% all animal summary analysis
% group by animal
close all
clear qqq wwFig rpc opc
[mgbTempTestsOnly,allDataTestsOnly,allDataCtlOnly]=byAnimalPercentCorrect(allDataTestsOnly,allDataCtlOnly,icDataTestsOnly,reinfcolor,optocolor);
close all
clear qqq wwFig rpc opc % now group by session, percent correct for Test
% [allDataTestsOnly,allDataCtlOnly]=bySessPercentCorrect(allDataTestsOnly,allDataCtlOnly,tempTestsOnly,mgbTempTestsOnly,reinfcolor,optocolor);
close all
 %% now do by animal for hit and FA for Test animals
% tempTestsOnly is the IC test animals, mgbTempTestsOnly is the MGB animals
mgbTempTestsOnly{1,39} = 'dprime MGB T Full, light on and light off';
mgbTempTestsOnly{1,40} = 'dprime MGB T Tone, light on and light off';
mgbTempTestsOnly{1,41} = 'dprime MGB T Choice, light on and light off';
allDataCtlOnly{1,39} = 'dprime MGB C Full, light on and light off';
allDataCtlOnly{1,40} = 'dprime MGB C Tone, light on and light off';
allDataCtlOnly{1,41} = 'dprime MGB C Choice, light on and light off';

close all
[mgbTempTestsOnly,allDataCtlOnly,allDataTestsOnly]=byAnimalHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,icDataTestsOnly,reinfcolor,optocolor);
close all
% bySessHFA(allDataTestsOnly,allDataCtlOnly,mgbTempTestsOnly,tempTestsOnly,reinfcolor,optocolor);
% close all
%%
% now do the anova....
[aovMGBbySess,statsMGBbySess,aovMGBbyAn,statsMGBbyAn,aovIC,statsIC]=anova2OptoPerAnimal(mgbTempTestsOnly,tempTestsOnly);
% close all
aovMGBbySess
aovMGBbyAn
% aovIC

%% compare differences in percent correct across each light off vs light on condition
%this is by session
% none of these analyses are significant 2025/01
% allDataTestsOnly{1,39}='MGB Difference PC Full Trial';
% allDataTestsOnly{1,40}='MGB Difference PC Tone';
% allDataTestsOnly{1,41}='MGB Difference PC Choice';
% allDataTestsOnly{1,42}='IC Difference PC Full Trial';
% allDataTestsOnly{1,43}='IC Difference PC Tone';
% allDataTestsOnly{1,44}='IC Difference PC Choice';
% % by session
% [allDataTestsOnly,mgbTempTestsOnly,tempTestsOnly]=diffPcOptoSess(allDataTestsOnly,mgbTempTestsOnly,tempTestsOnly,optocolor);
% % by animal
% diffPcOptoAn(allDataTestsOnly,optocolor);

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
mgbTempTestsOnly{1,42} = 'lick latency MGB T Full, hitLickRate; ohitLickRate; faLickRate; ofaLickRate';
mgbTempTestsOnly{1,43} = 'lick latency MGB T Tone, light on and light off';
mgbTempTestsOnly{1,44} = 'lick latency MGB T Choice, light on and light off';
mgbTempTestsOnly{1,45} = 'lick rate MGB T Full, light on and light off';
mgbTempTestsOnly{1,46} = 'lick rate MGB T Tone, light on and light off';
mgbTempTestsOnly{1,47} = 'lick rate MGB T Choice, light on and light off';

allDataCtlOnly{1,42} = 'lick latency MGB C Full, hitLickRate; ohitLickRate; faLickRate; ofaLickRate';
allDataCtlOnly{1,43} = 'lick latency MGB C Tone, light on and light off';
allDataCtlOnly{1,44} = 'lick latency MGB C Choice, light on and light off';
allDataCtlOnly{1,45} = 'lick rate MGB T Full, light on and light off';
allDataCtlOnly{1,46} = 'lick rate MGB T Tone, light on and light off';
allDataCtlOnly{1,47} = 'lick rate MGB T Choice, light on and light off';

lickLatOptoPerAnimal(mgbTempTestsOnly,allDataTestsOnly,allDataCtlOnly,reinfcolor,optocolor)

%% AC Expert Data
cd('O:\sjk\Behavior\AC9010Expert'); % this is celine's expert data where she inactivated on 90% of trials
expertACSummaryData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_day.mat'); 
expertACSummaryData=expertACSummaryData.summary_data;
expertACRawData=load('O:\sjk\Behavior\AC9010Expert\summary_90optoexpert_lickraster.mat');

plotACExpertData(expertACSummaryData);