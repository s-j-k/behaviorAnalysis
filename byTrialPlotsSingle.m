function byTrialPlotsSingle(allDataTestsOnly,allLicksTest,days,reinfcolor,optocolor)
% plot to make rolling behavior average
    counter=1;counter2=1;
    nbins=50;bins = linspace(-1,4,nbins);  
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
    
SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
    
lickIdx=3;
licks=allLicksTest{lickIdx,1};
% now get the relevant days
fullDays=days{lickIdx+1,4};
toneDays=days{lickIdx+1,5};
licksRFT=licks{fullDays(1,2)+1,1};
licksRFF=licks{fullDays(1,2)+1,2};
licksFT=licks{fullDays(1,2)+1,3};
licksFF=licks{fullDays(1,2)+1,4};
licksCFT=licks{fullDays(1,2)+1,7};
licksCFF=licks{fullDays(1,2)+1,8};

licksRTT=licks{toneDays(1,2)+1,1};
licksRTF=licks{toneDays(1,2)+1,2};
licksTT=licks{toneDays(1,2)+1,5};
licksTF=licks{toneDays(1,2)+1,6};
licksCTT=licks{toneDays(1,2)+1,7};
licksCTF=licks{toneDays(1,2)+1,8};

tempLickMat=NaN(300,60);
rFTLickMat=NaN(300,60);
rFFLickMat=NaN(300,60);
fTLickMat=NaN(300,60);
fFLickMat=NaN(300,60);
cFTLickMat=NaN(300,60);
cFFLickMat=NaN(300,60);
rRTLickMat=NaN(300,60);
rRFLickMat=NaN(300,60);
tTLickMat=NaN(300,60);
tFLickMat=NaN(300,60);
cTTLickMat=NaN(300,60);
cTFLickMat=NaN(300,60);
    % get consummatory licks
    for pp=1:2
        if pp==1
            dayIdx=fullDays(1,2);
        else
            dayIdx=toneDays(1,2);
        end
        nextIdx=1; %templicks is all of the lick data for that session;
        sessAllLicks=allLicksTest{lickIdx,2}{dayIdx};
        exampleSession=allDataTestsOnly{lickIdx+1,26}; 
        exampleTrials=find(exampleSession(:,1)==dayIdx);
        exampleTrials=exampleSession(exampleTrials,:);
        for tt=1:size(exampleTrials)-1
            nextIdxTemp=find(sessAllLicks>exampleTrials(tt,6));
            try
                nextIdx(tt+1)=nextIdxTemp(1);
            catch
                disp(tt)
            end
        end
        for yy=1:size(exampleTrials)-2
            tempLickMat(yy,1:length(sessAllLicks(nextIdx(yy):nextIdx(yy+1)-1)))=sessAllLicks(nextIdx(yy):nextIdx(yy+1)-1);
        end
        %now tempLickMat is the licks, for each trial, for the entire session
        %sort tempLickMat now by tone and context/condition
        if pp==1
            reinfTIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==1);
            reinfFIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==2);
            fullTIdx=find(exampleTrials(1:300,2)==1 & exampleTrials(1:300,3)==1);
            fullFIdx=find(exampleTrials(1:300,2)==1 & exampleTrials(1:300,3)==2);
            choiceTIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==1);
            choiceFIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==2);
            reinfFTLicks=tempLickMat(reinfTIdx,:);
            reinfFFLicks=tempLickMat(reinfFIdx,:);
            fullTLicks=tempLickMat(fullTIdx,:);
            fullFLicks=tempLickMat(fullFIdx,:);
            choiceFTLicks=tempLickMat(choiceTIdx,:);
            choiceFFLicks=tempLickMat(choiceFIdx,:);
        else
            reinfTIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==1);
            reinfFIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==2);
            toneTIdx=find(exampleTrials(1:300,2)==5 & exampleTrials(1:300,3)==1);
            toneFIdx=find(exampleTrials(1:300,2)==5 & exampleTrials(1:300,3)==2);
            choiceTIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==1);
            choiceFIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==2);

            reinfTTLicks=tempLickMat(reinfTIdx,:);
            reinfTFLicks=tempLickMat(reinfFIdx,:);
            toneTLicks=tempLickMat(toneTIdx,:);
            toneFLicks=tempLickMat(toneFIdx,:);
            choiceFTLicks=tempLickMat(choiceTIdx,:);
            choiceFFLicks=tempLickMat(choiceFIdx,:);
        end

        if pp==1
            reinfTIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==1);
            reinfFIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==2);
            fullTIdx=find(exampleTrials(1:300,2)==1 & exampleTrials(1:300,3)==1);
            fullFIdx=find(exampleTrials(1:300,2)==1 & exampleTrials(1:300,3)==2);
            choiceTIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==1);
            choiceFIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==2);
            for yy=1:size(exampleTrials)-2
                rFTLickMat(yy,1:length(licksRFT(nextIdx(yy):nextIdx(yy+1)-1)))=licksRFT(nextIdx(yy):nextIdx(yy+1)-1);
                rFFLickMat(yy,1:length(licksRFF(nextIdx(yy):nextIdx(yy+1)-1)))=licksRFF(nextIdx(yy):nextIdx(yy+1)-1);
                fTLickMat(yy,1:length(licksFT(nextIdx(yy):nextIdx(yy+1)-1)))=licksFT(nextIdx(yy):nextIdx(yy+1)-1);
                fFLickMat(yy,1:length(licksFF(nextIdx(yy):nextIdx(yy+1)-1)))=licksFF(nextIdx(yy):nextIdx(yy+1)-1);
                cFTLickMat(yy,1:length(licksCFT(nextIdx(yy):nextIdx(yy+1)-1)))=licksCFT(nextIdx(yy):nextIdx(yy+1)-1);
                cFFLickMat(yy,1:length(licksCFF(nextIdx(yy):nextIdx(yy+1)-1)))=licksCFF(nextIdx(yy):nextIdx(yy+1)-1);
            end
            rFTLickMat=rFTLickMat(reinfTIdx,:);
            rFFLickMat=rFFLickMat(reinfFIdx,:);
            fTLickMat=fTLickMat(fullTIdx,:);
            fFLickMat=fFLickMat(fullFIdx,:);
            cFTLickMat=cFTLickMat(choiceTIdx,:);
            cFFLickMat=cFFLickMat(choiceFIdx,:);
        else
            reinfTIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==1);
            reinfFIdx=find(exampleTrials(1:300,2)==2 & exampleTrials(1:300,3)==2);
            toneTIdx=find(exampleTrials(1:300,2)==5 & exampleTrials(1:300,3)==1);
            toneFIdx=find(exampleTrials(1:300,2)==5 & exampleTrials(1:300,3)==2);
            choiceTIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==1);
            choiceFIdx=find(exampleTrials(1:300,2)==6 & exampleTrials(1:300,3)==2);
            for yy=1:size(exampleTrials)-2
                rRTLickMat(yy,1:length(licksRTT(nextIdx(yy):nextIdx(yy+1)-1)))=licksRTT(nextIdx(yy):nextIdx(yy+1)-1);
                rRFLickMat(yy,1:length(licksRTF(nextIdx(yy):nextIdx(yy+1)-1)))=licksRTF(nextIdx(yy):nextIdx(yy+1)-1);
                tTLickMat(yy,1:length(licksTT(nextIdx(yy):nextIdx(yy+1)-1)))=licksTT(nextIdx(yy):nextIdx(yy+1)-1);
                tFLickMat(yy,1:length(licksTF(nextIdx(yy):nextIdx(yy+1)-1)))=licksTF(nextIdx(yy):nextIdx(yy+1)-1);
                cTTLickMat(yy,1:length(licksCTT(nextIdx(yy):nextIdx(yy+1)-1)))=licksCTT(nextIdx(yy):nextIdx(yy+1)-1);
                cTFLickMat(yy,1:length(licksCTF(nextIdx(yy):nextIdx(yy+1)-1)))=licksCTF(nextIdx(yy):nextIdx(yy+1)-1);
            end
            rRTLickMat=rRTLickMat(reinfTIdx,:);
            rRFLickMat=rRFLickMat(reinfFIdx,:);
            tTLickMat=tTLickMat(toneTIdx,:);
            tFLickMat=tFLickMat(toneFIdx,:);
            cTTLickMat=cTTLickMat(choiceTIdx,:);
            cTFLickMat=cTFLickMat(choiceFIdx,:);
        end
    end

    % now remove the probe block indices
    %'Condition','Trial indices','No probe indicies','Miss idx no probe'
    idxAllMice={'Reinf Tone T',allMiceRTTargetIdx;'Reinf Tone F',allMiceRTFoilIdx;...
        'Tone T',allMiceToneTargetIdx;'Tone F',allMiceToneFoilIdx;...
        'Choice Tone T',allMiceChoiceTTargetIdx;'Choice Tone F',allMiceChoiceTFoilIdx;...
        'Reinf Full T',allMiceRFTargetIdx;'Reinf Full F',allMiceRFFoilIdx;...
        'Full T',allMiceFullTargetIdx;'Full F',allMiceFullFoilIdx;...
        'Choice Full T',allMiceChoiceFTargetIdx;'Choice Full F',allMiceChoiceFFoilIdx};
    for yi=1:length(idxAllMice)
        for oo=1:size(idxAllMice{yi,2},1)
            probeIdx=find(idxAllMice{yi,2}(oo,:)>140);
            idxAllMice{yi,3}(oo,:)=[idxAllMice{yi,2}(oo,1:probeIdx(1)-1) ...
                (idxAllMice{yi,2}(oo,probeIdx(1):length(idxAllMice{yi,2}))-20)];
        end
    end
    
    idx=6;
    for rr=idx
        lickLatRFT=allMiceRFT{1,rr}(:,LICKL);
        lickLatRFF=allMiceRFF{1,rr}(:,LICKL);
        lickLatFullT=allMiceFullT{1,rr}(:,LICKL);
        lickLatFullF=allMiceFullF{1,rr}(:,LICKL);
        lickLatCFT=allMiceCFT{1,rr}(:,LICKL);
        lickLatCFF=allMiceCFF{1,rr}(:,LICKL);
        
        lickLatRFTHist=hist(lickLatRFT,bins)/length(lickLatRFT);
        lickLatRFFHist=hist(lickLatRFF,bins)/length(lickLatRFF);
        lickLatFullTHist=hist(lickLatFullT,bins)/length(lickLatFullT);
        lickLatFullFHist=hist(lickLatFullF,bins)/length(lickLatFullF);
        lickLatCFTHist=hist(lickLatCFT,bins)/length(lickLatCFT);
        lickLatCFFHist=hist(lickLatCFF,bins)/length(lickLatCFF);
        
        lickLatRTT=allMiceRTT{1,rr}(:,LICKL);
        lickLatRTF=allMiceRTF{1,rr}(:,LICKL);
        lickLatToneT=allMiceToneT{1,rr}(:,LICKL);
        lickLatToneF=allMiceToneF{1,rr}(:,LICKL);
        lickLatCT=allMiceCTT{1,rr}(:,LICKL);
        lickLatCF=allMiceCTF{1,rr}(:,LICKL);
        
        lickLatRTTHist=hist(lickLatRTT,bins)/length(lickLatRTT);
        lickLatRTFHist=hist(lickLatRTF,bins)/length(lickLatRTF);
        lickLatToneTHist=hist(lickLatToneT,bins)/length(lickLatToneT);
        lickLatToneFHist=hist(lickLatToneF,bins)/length(lickLatToneF);
        lickLatCTHist=hist(lickLatCT,bins)/length(lickLatCT);
        lickLatCFHist=hist(lickLatCF,bins)/length(lickLatCF);
    end
    
    lickLatRFTNoNan=lickLatRFT(~isnan(lickLatRFT));
    lickLatRFFNoNan=lickLatRFF(~isnan(lickLatRFF));    
    lickLatFullTNoNan=lickLatFullT(~isnan(lickLatFullT));
    lickLatFullFNoNan=lickLatFullF(~isnan(lickLatFullF));
    lickLatCFullTNoNan=lickLatCFT(~isnan(lickLatCFT));
    lickLatCFullFNoNan=lickLatCFF(~isnan(lickLatCFF));

    lickLatRTTNoNan=lickLatRTT(~isnan(lickLatRTT));
    lickLatTTNoNan=lickLatToneT(~isnan(lickLatToneT));
    lickLatRTFNoNan=lickLatRTF(~isnan(lickLatRTF));
    lickLatTFNoNan=lickLatToneF(~isnan(lickLatToneF));
    lickLatCTNoNan=lickLatCT(~isnan(lickLatCT));
    lickLatCFNoNan=lickLatCF(~isnan(lickLatCF));
    
    %add open circles aligned to 0 for the misses
    missIdx=find(allMiceRFTarget(idx,:)==0);
    missFIdx=find(allMiceFullTarget(idx,:)==0);
    missFCIdx=find(allMiceCFTarget(idx,:)==0);
    missRTIdx=find(allMiceRTTarget(idx,:)==0);
    missTIdx=find(allMiceToneTarget(idx,:)==0);
    missTCIdx=find(allMiceCTTarget(idx,:)==0);
    
    crIdx=find(allMiceRFFoil(idx,:)==0);
    crFIdx=find(allMiceFullFoil(idx,:)==0);
    crFCIdx=find(allMiceCFFoil(idx,:)==0);
    crRTIdx=find(allMiceRTFoil(idx,:)==0);
    crTIdx=find(allMiceToneFoil(idx,:)==0);
    crTCIdx=find(allMiceCTFoil(idx,:)==0);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% plots for scatter, PSTH %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idxShift=-0.5;
    sz=10;licksColor=[0.9 0.9 0.9];
    
    scatterPSTH=1;
    if scatterPSTH==1
        hitFullFig=figure;
        subplot(5,1,1);scatter(rFTLickMat(1:35,:),allMiceRFTargetIdx(idx,1:35)',sz,licksColor,'filled');
        hold on;scatter(lickLatRFT(1:35,:),allMiceRFTargetIdx(idx,1:35)',sz,reinfcolor,'filled');title('Light off, Hit');xlim([0 2]);
        ylabel('Trial');xlim([-0.5 3]);
        scatter(zeros(1,7)+idxShift,allMiceRFTargetIdx(idx,missIdx(1:7))',sz,[0.7 0 0]);

        subplot(5,1,2);    scatter(fTLickMat(1:17,:),allMiceFullTargetIdx(idx,1:17)',sz,licksColor,'filled');
        hold on;scatter(lickLatFullT(1:17,:),allMiceFullTargetIdx(idx,1:17)',sz,optocolor,'filled'); title('Full, Hit');xlim([0 2]);
        ylabel('Trial');xlim([-0.5 3]);
        scatter(zeros(1,8)+idxShift,allMiceFullTargetIdx(idx,missFIdx(1:8))',sz,'red');ylabel('Trial');

        subplot(5,1,3);scatter(cFTLickMat(1:17,:),allMiceChoiceFTargetIdx(idx,1:17)',sz,licksColor,'filled');
        hold on;scatter(lickLatCFT(1:17,:),allMiceChoiceFTargetIdx(idx,1:17)',sz,optocolor,'filled');title('Choice, Hit');xlim([0 2]);
        ylabel('Trial');xlabel('Time (s)');xlim([-0.5 3]);
        scatter(zeros(1,10)+idxShift,allMiceChoiceFTargetIdx(idx,missFCIdx(1:10))',sz,'red');ylabel('Trial');

        subplot(5,1,4);
        plot(bins,lickLatRFTHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatFullTHist,'Color',optocolor);xlabel('First lick latency (s)');
        ylabel('Proportion of S+ Trials');
        %     plot(1:length(lickLatRFTNoNan),lickLatRFTNoNan,'Color',reinfcolor);
        %     hold on;plot(1:length(lickLatFullTNoNan),lickLatFullTNoNan(1:length(lickLatFullTNoNan)),'Color','b');ylabel('First lick latency (s)');
        %     xlim([1 20]);
        box off;title('Full');
        subplot(5,1,5);ylabel('Proportion of S+ Trials');
        plot(bins,lickLatRFTHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatCFTHist,'Color',optocolor);
        %     plot(1:length(lickLatCFullTNoNan),lickLatRFTNoNan(1:length(lickLatCFullTNoNan)),'Color',reinfcolor);
        %     hold on;plot(1:length(lickLatCFullTNoNan),lickLatCFullTNoNan,'Color','b');xlim([1 20]);ylim([-0.5 4.5]);ylabel('First lick latency (s)');
        box off;title('Choice');
        xlabel('First lick latency (s)');
        hitFullFig.Position(3:4)=[250 875];
        saveas(gcf,['sk176D2_T_MGB_HitFull_OptoHist']);
        saveas(gcf,['sk176D2_T_MGB_HitFull_OptoHist.png']);    



        hitFig=figure;
        subplot(5,1,1);scatter(rRTLickMat(1:35,:),allMiceRTTargetIdx(idx,1:35)',sz,licksColor,'filled');
        hold on;scatter(lickLatRTT(1:35,:),allMiceRTTargetIdx(idx,1:35)',sz,reinfcolor,'filled');
        title('Light off, Hit');xlim([-0.5 3]);ylabel('Trial');
        %     scatter(zeros(1,length(missRTIdx)),allMiceRTTargetIdx(idx,missRTIdx)',sz,reinfcolor);
        subplot(5,1,2);scatter(tTLickMat(1:17,:),allMiceToneTargetIdx(idx,1:17)',sz,licksColor,'filled');
        hold on;scatter(lickLatToneT(1:17,:),allMiceToneTargetIdx(idx,1:17)',sz,optocolor,'filled'); title('Stimulus, Hit');xlim([0 2]);
        ylabel('Trial');scatter(zeros(1,5)+idxShift,allMiceToneTargetIdx(idx,missTIdx(1:5))',sz,'red');
        xlim([-0.5 3]);
        subplot(5,1,3);hold on;title('Stimulus, Hit');xlim([0 2]);
        scatter(cTTLickMat(1:17,:),allMiceChoiceTTargetIdx(idx,1:17)',sz,licksColor,'filled');
        scatter(lickLatCT(1:17,:),allMiceChoiceTTargetIdx(idx,1:17)',sz,optocolor,'filled'); 
        title('Choice, Hit');xlim([0 2]);
        ylabel('Trial');xlabel('Time (s)');hold on;
        scatter(zeros(1,8)+idxShift,allMiceChoiceTTargetIdx(idx,missTCIdx(1:8))',sz,'red');
        xlim([-0.5 3]);
        subplot(5,1,4);title('Stimulus');
        plot(bins,lickLatRTTHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatToneTHist,'Color',optocolor);xlabel('First lick latency (s)');
        ylabel('Proportion of S+ Trials');
        %     plot(1:length(lickLatTTNoNan),lickLatRTTNoNan(1:length(lickLatTTNoNan)),'Color',reinfcolor);
        %     hold on;plot(1:length(lickLatTTNoNan),lickLatTTNoNan,'Color','b');ylabel('First lick latency (s)');
        %     xlim([1 20]);box off;ylim([-0.5 4.5]);
        subplot(5,1,5);title('Choice');
        plot(bins,lickLatRTTHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatCTHist,'Color',optocolor);xlabel('First lick latency (s)');
        ylabel('Proportion of S+ Trials');
        %     plot(1:length(lickLatCTNoNan),lickLatTTNoNan(1:length(lickLatCTNoNan)),'Color',reinfcolor);
        %     hold on;plot(1:length(lickLatCTNoNan),lickLatCTNoNan,'Color','b');box off;
        %     xlabel('Hit Trial');ylabel('First lick latency (s)');xlim([1 20]);ylim([-0.5 4.5]);
        hitFig.Position(3:4)=[250 875];
        saveas(gcf,['sk176D2_T_MGB_Hit_OptoHist']);
        saveas(gcf,['sk176D2_T_MGB_Hit_OptoHist.png']);    

        FAFig=figure;
        subplot(5,1,1);hold on;scatter(rRFLickMat(1:35,:),allMiceRTFoilIdx(idx,1:35)',sz,licksColor,'filled');
        scatter(lickLatRTF(1:35,:),allMiceRTFoilIdx(idx,1:35)',sz,reinfcolor,'filled');title('Light off, FA');
        xlim([-0.5 3]);scatter(zeros(1,33)+idxShift,allMiceRTFoilIdx(idx,crIdx(1:33))',sz,'red');
        ylabel('Trial');
        subplot(5,1,2);scatter(tFLickMat(1:17,:),allMiceToneFoilIdx(idx,1:17)',sz,licksColor,'filled');
        hold on;scatter(lickLatToneF(1:17,:),allMiceToneFoilIdx(idx,1:17)',sz,optocolor,'filled'); title('Stimulus opto');
        xlim([-0.5 3]);scatter(zeros(1,6)+idxShift,allMiceToneFoilIdx(idx,crTIdx(1:6))',sz,'red');
        ylabel('Trial');
        subplot(5,1,3);scatter(cTFLickMat(1:17,:),allMiceChoiceTFoilIdx(idx,1:17)',sz,licksColor,'filled');
        hold on;scatter(lickLatCF(1:17,:),allMiceChoiceTFoilIdx(idx,1:17)',sz,optocolor,'filled');title('Choice opto');
        xlim([-0.5 3]);scatter(zeros(1,12)+idxShift,allMiceChoiceTFoilIdx(idx,crTCIdx(1:12))',sz,'red');
        ylabel('Trial');xlabel('Time (s)');
        subplot(5,1,4);
        plot(bins,lickLatRTFHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatFullFHist,'Color',optocolor);xlabel('First lick latency (s)');
        ylabel('Proportion of S- Trials');
        %     shadedErrorBar(1:length(lickLatRTFNoNan),lickLatRTFNoNan(1:length(lickLatRTFNoNan)),'Color',reinfcolor);
        %     hold on;shadedErrorBar(1:length(lickLatTFNoNan),lickLatTFNoNan,'Color','b');ylabel('First lick latency (s)');
        %     xlim([1 20]);box off;
        subplot(5,1,5);
        plot(bins,lickLatRTFHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatCFHist,'Color',optocolor);xlabel('First lick latency (s)');
        ylabel('Proportion of S- Trials');
        %     plot(1:length(lickLatCFNoNan),lickLatRFTNoNan(1:length(lickLatCFNoNan)),'Color',reinfcolor);
        %     hold on;plot(1:length(lickLatCFNoNan),lickLatCFNoNan,'Color','b');box off;
        %     xlabel('FA Trial');ylabel('First lick latency (s)');xlim([1 20]);
        FAFig.Position(3:4)=[250 875];
        saveas(gcf,['sk176D2_T_MGB_FA_OptoHist']);
        saveas(gcf,['sk176D2_T_MGB_FA_OptoHist.png']);    


        FAFullFig=figure;
        subplot(5,1,1);hold on;scatter(rFFLickMat(1:35,:),allMiceRFFoilIdx(idx,1:35)',sz,licksColor,'filled');
        scatter(lickLatRFF(1:35,:),allMiceRFFoilIdx(idx,1:35)',sz,reinfcolor,'filled');title('Light off, FA');
        ylabel('Trial');scatter(zeros(1,33)+idxShift,allMiceRFFoilIdx(idx,crIdx(1:33))',sz,'red');xlim([-0.5 3])
        subplot(5,1,2);hold on;scatter(fFLickMat(1:17,:),allMiceFullFoilIdx(idx,1:17)',sz,licksColor,'filled');
        scatter(lickLatFullF(1:17),allMiceFullFoilIdx(idx,1:17)',sz,optocolor,'filled'); title('Full');xlim([0 2]);
        scatter(zeros(1,13)+idxShift,allMiceFullFoilIdx(idx,crFIdx(1:13))',sz,'red');xlim([-0.5 3])
        ylabel('Trial');
        subplot(5,1,3);hold on;scatter(cFFLickMat(1:17,:),allMiceChoiceFFoilIdx(idx,1:17)',sz,licksColor,'filled');
        scatter(lickLatCFF(1:17,:),allMiceChoiceFFoilIdx(idx,1:17)',sz,optocolor,'filled');title('Choice');
        scatter(zeros(1,11)+idxShift,allMiceChoiceFFoilIdx(idx,crFCIdx(1:11))',sz,'red');xlim([-0.5 3])
        ylabel('Trial');xlabel('Time (s)');
        subplot(5,1,4);
        plot(bins,lickLatRFFHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatFullFHist,'Color',optocolor);xlabel('First lick latency (s)');
        ylabel('Proportion of S- Trials');
        %     plot(1:length(lickLatRFFNoNan),lickLatRFFNoNan,'Color',reinfcolor);
        %     hold on;plot(1:length(lickLatFullFNoNan),lickLatFullFNoNan,'b');ylabel('First lick latency (s)');
        %     xlim([1 10]);box off;
        subplot(5,1,5);
        plot(bins,lickLatRFFHist,'Color',reinfcolor); hold on;
        plot(bins,lickLatCFFHist,'Color',optocolor);xlabel('First lick latency (s)');
        ylabel('Proportion of S- Trials');
        %     plot(1:length(lickLatRFFNoNan),lickLatRFFNoNan,'Color',reinfcolor);
        %     hold on;plot(1:length(lickLatCFullFNoNan),lickLatCFullFNoNan,'b');box off;
        %     xlabel('FA Trial');ylabel('First lick latency (s)');xlim([1 10]);
        FAFullFig.Position(3:4)=[250 875];
        saveas(gcf,['sk176D2_T_MGB_Example_AnimalFAFull_OptoHist']);
        saveas(gcf,['sk176D2_T_MGB_Example_AnimalFAFull_OptoHist.png']);    

    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% plots for rolling average %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rollingb=0;
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

        %%%%%scatter plot
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
    
    
    
end