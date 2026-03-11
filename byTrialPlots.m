function byTrialPlots(allDataTestsOnly,days,reinfcolor,optocolor)
% plot to make rolling behavior average
    counter=1;counter2=1;
    for qq=2:size(allDataTestsOnly,1)
        exampleSession=allDataTestsOnly{qq,26}; 
        mgbDays=days{qq,5}; 
        for ww=1:2
            sessionId = mgbDays(ww); % tone day 1
            exampleTrials=find(exampleSession(:,1)==sessionId);
            exampleTrials=exampleSession(exampleTrials,:);

            reinfDataD1TargetIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==1);
            reinfDataD1FoilIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==2);
            reinfDataD1Target=exampleTrials(reinfDataD1TargetIdx,:);
            reinfDataD1Foil=exampleTrials(reinfDataD1FoilIdx,:);
            % now relabel misses as zeroes instead of twos
            % TXT = 2; TONE = 3; OUTCOME = 4; 
            missIdx=find(reinfDataD1Target(:,4)==2);
            reinfDataD1Target(missIdx,4)=0;
            %relabel CR and FA as 0 and 1
            crIdx=find(reinfDataD1Foil(:,4)==4);
            reinfDataD1Foil(crIdx,4)=0;
            faIdx=find(reinfDataD1Foil(:,4)==3);
            reinfDataD1Foil(faIdx,4)=1;

            %Tone
            toneDataD1TargetIdx=find(exampleTrials(1:300,2)==5 & exampleTrials(1:300,3)==1);
            toneDataD1FoilIdx=find(exampleTrials(1:300,2)==5 & exampleTrials(1:300,3)==2);
            toneDataD1Target=exampleTrials(toneDataD1TargetIdx,:);
            toneDataD1Foil=exampleTrials(toneDataD1FoilIdx,:);
            % now relabel false alarms as zeroes instead of twos
            missIdx=find(toneDataD1Target(:,4)==2);
            toneDataD1Target(missIdx,4)=0;
            %relabel CR and FA as 0 and 1
            crIdx=find(toneDataD1Foil(:,4)==4);
            toneDataD1Foil(crIdx,4)=0;
            faIdx=find(toneDataD1Foil(:,4)==3);
            toneDataD1Foil(faIdx,4)=1;

            choiceDataD1TargetIdx=find(exampleTrials(:,2)==6 & exampleTrials(:,3)==1);
            choiceDataD1FoilIdx=find(exampleTrials(:,2)==6 & exampleTrials(:,3)==2);
            choiceDataD1Target=exampleTrials(choiceDataD1TargetIdx,:);
            choiceDataD1Foil=exampleTrials(choiceDataD1FoilIdx,:);
            % now relabel false alarms as zeroes instead of twos
            missIdx=find(choiceDataD1Target(:,4)==2);
            choiceDataD1Target(missIdx,4)=0;
            %relabel CR and FA as 0 and 1
            crIdx=find(choiceDataD1Foil(:,4)==4);
            choiceDataD1Foil(crIdx,4)=0;
            faIdx=find(choiceDataD1Foil(:,4)==3);
            choiceDataD1Foil(faIdx,4)=1;
            
            allMiceRTTargetIdx(counter,:)=reinfDataD1TargetIdx;
            allMiceRTFoilIdx(counter,:)=reinfDataD1FoilIdx;
            allMiceToneTargetIdx(counter,:)=toneDataD1TargetIdx;
            allMiceToneFoilIdx(counter,:)=toneDataD1FoilIdx;
            allMiceChoiceTTargetIdx(counter,:)=choiceDataD1TargetIdx;
            allMiceChoiceTFoilIdx(counter,:)=choiceDataD1FoilIdx;
            
            allMiceRTTarget(counter,:)=reinfDataD1Target(1:70,4);
            allMiceRTFoil(counter,:)=reinfDataD1Foil(1:70,4);
            allMiceToneTarget(counter,:)=toneDataD1Target(1:35,4);
            allMiceToneFoil(counter,:)=toneDataD1Foil(1:35,4);
            allMiceCTTarget(counter,:)=choiceDataD1Target(1:35,4);
            allMiceCTFoil(counter,:)=choiceDataD1Foil(1:35,4);
            
            allMiceRTT{counter}=reinfDataD1Target;
            allMiceRTF{counter}=reinfDataD1Foil;
            allMiceToneT{counter}=toneDataD1Target;
            allMiceToneF{counter}=toneDataD1Foil;
            allMiceCTT{counter}=choiceDataD1Target;
            allMiceCTF{counter}=choiceDataD1Foil;
            counter=counter+1;
        end

        % full trial data
        mgbDays=days{qq,4}; % full trial days
        clear reinfDataD1 reinfDataD1Target reinfDataD1Foil
        for ee=1:2
            sessionId = mgbDays(ee); 
            exampleTrials=find(exampleSession(:,1)==sessionId);
            exampleTrials=exampleSession(exampleTrials,:);
            reinfDataD1TargetIdx=find(exampleTrials(:,2)==2 & exampleTrials(:,3)==1);
            reinfDataD1FoilIdx=find(exampleTrials(:,2)==2 & exampleTrials(:,3)==2);
            reinfDataD1Target=exampleTrials(reinfDataD1TargetIdx,:);
            reinfDataD1Foil=exampleTrials(reinfDataD1FoilIdx,:);
            % now relabel misses as zeroes instead of twos
            missIdx=find(reinfDataD1Target(:,4)==2);
            reinfDataD1Target(missIdx,4)=0;
            %relabel CR and FA as 0 and 1
            crIdx=find(reinfDataD1Foil(:,4)==4);
            reinfDataD1Foil(crIdx,4)=0;
            faIdx=find(reinfDataD1Foil(:,4)==3);
            reinfDataD1Foil(faIdx,4)=1;


            %Full trial
            fullDataD1TargetIdx=find(exampleTrials(1:300,2)==1 & exampleTrials(1:300,3)==1);
            fullDataD1FoilIdx=find(exampleTrials(1:300,2)==1 & exampleTrials(1:300,3)==2);
            fullDataD1Target=exampleTrials(fullDataD1TargetIdx,:);
            fullDataD1Foil=exampleTrials(fullDataD1FoilIdx,:);
            missIdx=find(fullDataD1Target(:,4)==2);
            fullDataD1Target(missIdx,4)=0;
            %relabel CR and FA as 0 and 1
            crIdx=find(fullDataD1Foil(:,4)==4);
            fullDataD1Foil(crIdx,4)=0;
            faIdx=find(fullDataD1Foil(:,4)==3);
            fullDataD1Foil(faIdx,4)=1;
            
            
            choiceDataD1TargetIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==1);
            choiceDataD1FoilIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==2);
            choiceDataD1Target=exampleTrials(choiceDataD1TargetIdx,:);
            choiceDataD1Foil=exampleTrials(choiceDataD1FoilIdx,:);
            missIdx=find(choiceDataD1Target(:,4)==2);
            choiceDataD1Target(missIdx,4)=0;
            %relabel CR and FA as 0 and 1
            crIdx=find(choiceDataD1Foil(:,4)==4);
            choiceDataD1Foil(crIdx,4)=0;
            faIdx=find(choiceDataD1Foil(:,4)==3);
            choiceDataD1Foil(faIdx,4)=1;
            
            allMiceRFTargetIdx(counter2,:)=reinfDataD1TargetIdx;
            allMiceRFFoilIdx(counter2,:)=reinfDataD1FoilIdx;
            allMiceFullTargetIdx(counter2,:)=toneDataD1TargetIdx;
            allMiceFullFoilIdx(counter2,:)=toneDataD1FoilIdx;
            allMiceChoiceFTargetIdx(counter2,:)=choiceDataD1TargetIdx;
            allMiceChoiceFFoilIdx(counter2,:)=choiceDataD1FoilIdx;
            
            allMiceRFTarget(counter2,:)=reinfDataD1Target(1:70,4);
            allMiceRFFoil(counter2,:)=reinfDataD1Foil(1:70,4);
            allMiceFullTarget(counter2,:)=fullDataD1Target(1:35,4);
            allMiceFullFoil(counter2,:)=fullDataD1Foil(1:35,4);
            allMiceCFTarget(counter2,:)=choiceDataD1Target(1:35,4);
            allMiceCFFoil(counter2,:)=choiceDataD1Foil(1:35,4);
            
            allMiceRFT{counter2}=reinfDataD1Target;
            allMiceRFF{counter2}=reinfDataD1Foil;
            allMiceFullT{counter2}=fullDataD1Target;
            allMiceFullF{counter2}=fullDataD1Foil;
            allMiceCFT{counter2}=choiceDataD1Target;
            allMiceCFF{counter2}=choiceDataD1Foil;
            
            counter2=counter2+1;
        end    
    end
    
%lick latency by trial 
SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
    
idx=6;
for rr=idx
    lickLatRFT=allMiceRFT{1,rr}(:,LICKL);
    lickLatRFF=allMiceRFF{1,rr}(:,LICKL);
    lickLatFullT=allMiceFullT{1,rr}(:,LICKL);
    lickLatFullF=allMiceFullF{1,rr}(:,LICKL);
    lickLatCFT=allMiceCFT{1,rr}(:,LICKL);
    lickLatCFF=allMiceCFF{1,rr}(:,LICKL);

    lickLatRTT=allMiceRTT{1,rr}(:,LICKL);
    lickLatRTF=allMiceRTF{1,rr}(:,LICKL);
    lickLatToneT=allMiceToneT{1,rr}(:,LICKL);
    lickLatToneF=allMiceToneF{1,rr}(:,LICKL);
    lickLatCT=allMiceCTT{1,rr}(:,LICKL);
    lickLatCF=allMiceCTF{1,rr}(:,LICKL);
end
    
lickLatRFTNoNan=lickLatRFT(~isnan(lickLatRFT));
lickLatRFFNoNan=lickLatRFF(~isnan(lickLatRFF));    
lickLatFullTNoNan=lickLatFullT(~isnan(lickLatFullT));
lickLatFullTFNoNan=lickLatFullF(~isnan(lickLatFullF));
lickLatCFTNoNan=lickLatCFT(~isnan(lickLatCFT));
lickLatCFFNoNan=lickLatCFF(~isnan(lickLatCFF));

lickLatRTTNoNan=lickLatRTT(~isnan(lickLatRTT));
lickLatTTNoNan=lickLatToneT(~isnan(lickLatToneT));
lickLatRTFNoNan=lickLatRTF(~isnan(lickLatRTF));
lickLatTFNoNan=lickLatToneF(~isnan(lickLatToneF));
lickLatCTNoNan=lickLatCT(~isnan(lickLatCT));
lickLatCFNoNan=lickLatCF(~isnan(lickLatCF));
    
rollingb=1;
    if rollingb==1
        % averaged across all mice plots
        windowsize=10;
        figure; subplot(3,1,1); title('Light off'); hold on
        stdTTarget=std(movmean(mean(allMiceRTTarget,1),windowsize));
        stdTFoil=std(movmean(mean(allMiceRTFoil,1),windowsize));
        shadedErrorBar([],movmean(mean(allMiceRTTarget,1),windowsize),std(movmean(mean(allMiceRTTarget,1),windowsize)),'g',0);hold on
        axis tight; ylim([0 1]);
        shadedErrorBar([],movmean(mean(allMiceRTFoil,1),windowsize),std(movmean(mean(allMiceRTFoil,1),windowsize)), 'r',0); legend('','hit','','fa');
        ylabel(['Rate, rolling average of ' num2str(windowsize) ' trials']);

        stdToneTarget=std(movmean(mean(allMiceToneTarget,1),windowsize));
        stdToneFoil=std(movmean(mean(allMiceToneFoil,1),windowsize));
        subplot(3,1,2);title('Stimulus');hold on
        shadedErrorBar([],movmean(mean(allMiceToneTarget),windowsize),std(movmean(mean(allMiceToneTarget,1),windowsize)),'g',0);hold on
        shadedErrorBar([],movmean(mean(allMiceToneFoil),windowsize),std(movmean(mean(allMiceToneFoil,1),windowsize)),'r',0);ylim([0 1]);
        axis tight; ylim([0 1]);

        stdChoiceTTarget=std(movmean(mean(allMiceCTTarget,1),windowsize));
        stdChoiceTFoil=std(movmean(mean(allMiceCTFoil,1),windowsize));
        subplot(3,1,3); title('Choice');hold on
        shadedErrorBar([],movmean(mean(allMiceCTTarget),windowsize),std(movmean(mean(allMiceCTTarget,1),windowsize)),'g',0);hold on
        shadedErrorBar([],movmean(mean(allMiceCTFoil),windowsize),std(movmean(mean(allMiceToneFoil,1),windowsize)),'r',0);
        axis tight; ylim([0 1]);
        xlabel('Trials');

        stdFTarget=std(movmean(mean(allMiceRFTarget,1),windowsize));
        stdFFoil=std(movmean(mean(allMiceRFFoil,1),windowsize));
        figure; subplot(3,1,1); title('Light off'); hold on
        shadedErrorBar([],movmean(mean(allMiceRFTarget),windowsize),std(movmean(mean(allMiceRFTarget,1),windowsize)),'g',0);hold on
        shadedErrorBar([],movmean(mean(allMiceRFFoil),windowsize),std(movmean(mean(allMiceRFFoil,1),windowsize)),'r',0); legend('','hit','','fa');
        axis tight; ylim([0 1]);
        ylabel(['Rate, rolling average of ' num2str(windowsize) ' trials']);

        stdFullTarget=std(movmean(mean(allMiceFullTarget,1),windowsize));
        stdFullFoil=std(movmean(mean(allMiceFullFoil,1),windowsize));
        subplot(3,1,2);title('Full trial');hold on
        shadedErrorBar([],movmean(mean(allMiceFullTarget),windowsize),std(movmean(mean(allMiceFullTarget,1),windowsize)),'g',0);hold on
        shadedErrorBar([],movmean(mean(allMiceFullFoil),windowsize),std(movmean(mean(allMiceFullFoil,1),windowsize)),'r',0); ylim([0 1]);
        axis tight; ylim([0 1]);

        stdChoiceTTarget=std(movmean(mean(allMiceCFTarget,1),windowsize));
        stdChoiceTFoil=std(movmean(mean(allMiceCFFoil,1),windowsize));   
        subplot(3,1,3); title('Choice');hold on
        shadedErrorBar([],movmean(mean(allMiceCFTarget),windowsize),std(movmean(mean(allMiceCFTarget,1),windowsize)),'g',0);hold on
        shadedErrorBar([],movmean(mean(allMiceCFFoil),windowsize),std(movmean(mean(allMiceCFFoil,1),windowsize)),'r',0);
        axis tight; ylim([0 1]);
        xlabel('Trials');
        ylim([0 1]);

        %scatter plot
        targetIdx = movmean(mean(allMiceRTTargetIdx,1),windowsize)';
        targetData=movmean(mean(allMiceRTTarget,1),windowsize)';
        foilIdx=movmean(mean(allMiceRTFoilIdx,1),windowsize)';
        foilData=movmean(mean(allMiceRTFoil,1),windowsize)';
        windowsize=10;
        figure; subplot(3,1,1); title('Light off'); hold on
        scatter(targetIdx,targetData);hold on;ylim([0 1]);
        scatter(foilIdx,foilData); legend('hit','fa');
        axis tight; ylim([0 1]);
        ylabel(['Rate, rolling average of ' num2str(windowsize) ' trials']);

        toneIdx = movmean(mean(allMiceToneTargetIdx,1),windowsize)';
        toneData=movmean(mean(allMiceToneTarget,1),windowsize)';
        toneFoilIdx=movmean(mean(allMiceToneFoilIdx,1),windowsize)';
        toneFoilData=movmean(mean(allMiceToneFoil,1),windowsize)';
        subplot(3,1,2);title('Stimulus');hold on
        scatter(toneIdx,toneData);hold on;ylim([0 1]);
        scatter(toneFoilIdx,toneFoilData); 
        axis tight; ylim([0 1]);

        choiceIdx=movmean(mean(allMiceChoiceTTargetIdx,1),windowsize)';
        choiceData=movmean(mean(allMiceCTTarget,1),windowsize)';
        choiceFoilIdx=movmean(mean(allMiceChoiceTFoilIdx,1),windowsize)';
        choiceFoilData=movmean(mean(allMiceCTFoil,1),windowsize)';
        subplot(3,1,3); title('Choice');hold on
        scatter(choiceIdx,choiceData);hold on;ylim([0 1]);
        scatter(choiceFoilIdx,choiceFoilData); 
        axis tight; ylim([0 1]);xlabel('Trials');

        rfTargetIdx=movmean(mean(allMiceRFTargetIdx,1),windowsize)';
        rfTargetData=movmean(mean(allMiceRFTarget,1),windowsize)';
        rfFoilIdx=movmean(mean(allMiceRFFoilIdx,1),windowsize)';
        rfFoilData=movmean(mean(allMiceRFFoil,1),windowsize)';
        figure; subplot(3,1,1); title('Light off'); hold on
        scatter(rfTargetIdx,rfTargetData);hold on;ylim([0 1]);
        scatter(rfFoilIdx,rfFoilData); legend('hit','fa');
        axis tight; ylim([0 1]);
        ylabel(['Rate, rolling average of ' num2str(windowsize) ' trials']);

        fullTargetIdx=movmean(mean(allMiceFullTargetIdx,1),windowsize)';
        fullTargetData=movmean(mean(allMiceFullTarget,1),windowsize)';
        fullFoilIdx=movmean(mean(allMiceFullFoilIdx,1),windowsize)';
        fullFoilData=movmean(mean(allMiceFullFoil,1),windowsize)';
        subplot(3,1,2);title('Full trial');hold on
        scatter(fullTargetIdx,fullTargetData);hold on;ylim([0 1]);
        scatter(fullFoilIdx,fullFoilData); 
        axis tight; ylim([0 1]);

        cfTargetIdx=movmean(mean(allMiceChoiceFTargetIdx,1),windowsize)';
        cfTargetData=movmean(mean(allMiceCFTarget,1),windowsize)';
        cfFoilIdx=movmean(mean(allMiceChoiceFFoilIdx,1),windowsize)';
        cfFoilData=movmean(mean(allMiceCFFoil,1),windowsize)';
        subplot(3,1,3); title('Choice');hold on
        scatter(cfTargetIdx,cfTargetData);hold on;ylim([0 1]);
        scatter(cfFoilIdx,cfFoilData); 
        axis tight; ylim([0 1]);
        xlabel('Trials');

    else
    end
    
    %add open circles aligned to 0 for the misses
    missIdx=find(allMiceRFTarget(idx,:)==0);
    missFIdx=find(allMiceFullTarget(idx,:)==0);
    missFCIdx=find(allMiceCFTarget(idx,:)==0);
    missRTIdx=find(allMiceRTTarget(idx,:)==0);
    missTIdx=find(allMiceToneTarget(idx,:)==0);
    missTCIdx=find(allMiceCTTarget(idx,:)==0);
    % do this also for teh CRs
    % add consummatory licks
    sz=20;
    hitFullFig=figure;
    subplot(5,1,1);scatter(lickLatRFT,allMiceRFTargetIdx(idx,:)',sz,reinfcolor,'filled');title('Light off, Hit');xlim([0 2]);
    hold on;ylabel('Trial');
    scatter(zeros(1,length(missIdx)),allMiceRFTargetIdx(idx,missIdx)',sz,reinfcolor);
    
    subplot(5,1,2);scatter(lickLatFullT,allMiceFullTargetIdx(idx,:)',sz,optocolor,'filled'); title('Full opto');xlim([0 2]);
    hold on;ylabel('Trial');ylim([-0.5 4.5]);
    scatter(zeros(1,length(missFIdx)),allMiceFullTargetIdx(idx,missFIdx)',sz,optocolor);ylabel('Trial');
    
    subplot(5,1,3);scatter(lickLatCT,allMiceChoiceFTargetIdx(idx,:)',sz,optocolor,'filled');title('Choice opto');xlim([0 2]);
    ylabel('Trial');xlabel('Time (s)');ylim([-0.5 4.5]);
    hold on;scatter(zeros(1,length(missFCIdx)),allMiceChoiceFTargetIdx(idx,missFCIdx)',sz,optocolor);ylabel('Trial');
    
    subplot(5,1,4);shadedErrorBar(1:length(lickLatRFTNoNan),lickLatRFTNoNan,std(lickLatRFTNoNan));
    hold on;shadedErrorBar(1:length(lickLatFullTNoNan),lickLatFullTNoNan(1:length(lickLatFullTNoNan)),std(lickLatFullTNoNan),'b');ylabel('First lick latency (s)');
    xlim([1 20]);box off;
    subplot(5,1,5);shadedErrorBar(1:length(lickLatCFTNoNan),lickLatRFTNoNan(1:length(lickLatCFTNoNan)),std(lickLatRFTNoNan));
    hold on;shadedErrorBar(1:length(lickLatCFTNoNan),lickLatCFTNoNan,std(lickLatCFTNoNan),'b');box off;
    xlabel('Hit Trial');ylabel('First lick latency (s)');xlim([1 20]);
    hitFullFig.Position(3:4)=[250 875];
    saveas(gcf,['sk176D2_T_MGB_Example_AnimalHitFull_Opto']);
    saveas(gcf,['sk176D2_T_MGB_Example_AnimalHitFull_Opto.png']);    
    
    hitFig=figure;
    subplot(5,1,1);scatter(lickLatRFT,allMiceRTTargetIdx(idx,:)',sz,reinfcolor,'filled');title('Light off, Hit');xlim([0 2]);
    ylabel('Trial');hold on;
    scatter(zeros(1,length(missRTIdx)),allMiceRTTargetIdx(idx,missRTIdx)',sz,reinfcolor);
    subplot(5,1,2);scatter(lickLatToneT,allMiceToneTargetIdx(idx,:)',sz,optocolor,'filled'); title('Stimulus opto');xlim([0 2]);
    ylabel('Trial');hold on;scatter(zeros(1,length(missTIdx)),allMiceToneTargetIdx(idx,missTIdx)',sz,optocolor);
    
    subplot(5,1,3);scatter(lickLatCT,allMiceChoiceTTargetIdx(idx,:)',sz,optocolor,'filled');title('Choice opto');xlim([0 2]);
    ylabel('Trial');xlabel('Time (s)');hold on;
    scatter(zeros(1,length(missTCIdx)),allMiceChoiceTTargetIdx(idx,missTCIdx)',sz,optocolor);
    subplot(5,1,4);shadedErrorBar(1:length(lickLatTTNoNan),lickLatRTTNoNan(1:length(lickLatTTNoNan)),std(lickLatRTTNoNan));
    hold on;shadedErrorBar(1:length(lickLatTTNoNan),lickLatTTNoNan,std(lickLatTTNoNan),'b');ylabel('First lick latency (s)');
    xlim([1 20]);box off;ylim([-0.5 4.5]);
    subplot(5,1,5);shadedErrorBar(1:length(lickLatCTNoNan),lickLatTTNoNan(1:length(lickLatCTNoNan)),std(lickLatTTNoNan));
    hold on;shadedErrorBar(1:length(lickLatCTNoNan),lickLatCTNoNan,std(lickLatCTNoNan),'b');box off;
    xlabel('Hit Trial');ylabel('First lick latency (s)');xlim([1 20]);ylim([-0.5 4.5]);
    hitFig.Position(3:4)=[250 875];
    saveas(gcf,['sk176D2_T_MGB_Example_AnimalHit_Opto']);
    saveas(gcf,['sk176D2_T_MGB_Example_AnimalHit_Opto.png']);    
    
    FAFig=figure;
    subplot(5,1,1);scatter(lickLatRTF,allMiceRTFoilIdx(idx,:)',sz,reinfcolor,'filled');title('Light off, FA');xlim([0 2]);ylim([1 300]);
    ylabel('Trial');
    subplot(5,1,2);scatter(lickLatToneF,allMiceToneFoilIdx(idx,:)',sz,optocolor,'filled'); title('Stimulus opto');xlim([0 2]);ylim([1 300]);
    ylabel('Trial');
    subplot(5,1,3);scatter(lickLatCF,allMiceChoiceTFoilIdx(idx,:)',sz,optocolor,'filled');title('Choice opto');xlim([0 2]);ylim([1 300]);
    ylabel('Trial');xlabel('Time (s)');
    subplot(5,1,4);shadedErrorBar(1:length(lickLatRTFNoNan),lickLatRTFNoNan(1:length(lickLatRTFNoNan)),std(lickLatRTFNoNan));
    hold on;shadedErrorBar(1:length(lickLatTFNoNan),lickLatTFNoNan,std(lickLatRTFNoNan),'b');ylabel('First lick latency (s)');
    xlim([1 20]);box off;
    subplot(5,1,5);shadedErrorBar(1:length(lickLatCFNoNan),lickLatRFTNoNan(1:length(lickLatCFNoNan)),std(lickLatRFTNoNan));
    hold on;shadedErrorBar(1:length(lickLatCFNoNan),lickLatCFNoNan,std(lickLatCFNoNan),'b');box off;
    xlabel('FA Trial');ylabel('First lick latency (s)');xlim([1 20]);
    FAFig.Position(3:4)=[250 875];
    saveas(gcf,['sk176D2_T_MGB_Example_AnimalFA_Opto']);
    saveas(gcf,['sk176D2_T_MGB_Example_AnimalFA_Opto.png']);    
    
    
end
