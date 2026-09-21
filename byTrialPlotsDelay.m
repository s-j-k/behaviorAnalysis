function byTrialPlotsDelay(allDataTestsOnly,optomeanMat,allLicksTest,reinfcolor,optocolor)
CTXT=2;

% Build the delay matrices in mouse/session order. Each row in a given
% delay matrix represents one unique session for one mouse.
    delayBySession={...
        [0,5,1],[0,2,1],[0,2,1];...
        [3,4,0],[3,4,0],[3,4,0];...
        [5,1,0],[5,1,0],[5,2,0];...
        [5,0,2],[5,2,0],[5,1,0];...
        [3,4,0],[3,4,0],[3,0,4]};
    contextOrder=[1,5,6];
    delayCounter=ones(1,5);

    delayRTargetIdx=cell(1,5); delayRFoilIdx=cell(1,5);
    delayTargetIdx=cell(1,5);  delayFoilIdx=cell(1,5);
    delayRTarget=cell(1,5);    delayRFoil=cell(1,5);
    delayTarget=cell(1,5);     delayFoil=cell(1,5);
    delayRT=cell(1,5);         delayRF=cell(1,5);
    delayT=cell(1,5);          delayF=cell(1,5);
    delaySessionKey=cell(1,5);

    for nbsubj=2:size(optomeanMat,1)
        exampleSession=optomeanMat{nbsubj,17};
        mouseIdx=nbsubj-1;

        for sessionId=1:size(delayBySession,1)
            exampleTrials=exampleSession(exampleSession(:,1)==sessionId,:);

            reinfTargetIdx=find(exampleTrials(:,CTXT)==2 & exampleTrials(:,3)==1);
            reinfFoilIdx=find(exampleTrials(:,CTXT)==2 & exampleTrials(:,3)==2);
            reinfTarget=exampleTrials(reinfTargetIdx,:);
            reinfFoil=exampleTrials(reinfFoilIdx,:);
            reinfTarget(reinfTarget(:,4)==2,4)=0;
            reinfFoil(reinfFoil(:,4)==4,4)=0;
            reinfFoil(reinfFoil(:,4)==3,4)=1;

            sessionDelays=delayBySession{sessionId,mouseIdx};
            for conditionIdx=1:numel(contextOrder)
                delayNumber=sessionDelays(conditionIdx);
                if delayNumber==0
                    continue
                end

                conditionContext=contextOrder(conditionIdx);
                targetIdx=find(exampleTrials(1:300,CTXT)==conditionContext & ...
                    exampleTrials(1:300,3)==1);
                foilIdx=find(exampleTrials(1:300,CTXT)==conditionContext & ...
                    exampleTrials(1:300,3)==2);
                target=exampleTrials(targetIdx,:);
                foil=exampleTrials(foilIdx,:);
                target(target(:,4)==2,4)=0;
                foil(foil(:,4)==4,4)=0;
                foil(foil(:,4)==3,4)=1;

                row=delayCounter(delayNumber);
                delaySessionKey{delayNumber}(row,:)=[mouseIdx,sessionId];
                delayRTargetIdx{delayNumber}(row,:)=reinfTargetIdx;
                delayRFoilIdx{delayNumber}(row,:)=reinfFoilIdx;
                delayTargetIdx{delayNumber}(row,:)=targetIdx;
                delayFoilIdx{delayNumber}(row,:)=foilIdx;
                delayRTarget{delayNumber}(row,:)=reinfTarget(1:70,4);
                delayRFoil{delayNumber}(row,:)=reinfFoil(1:70,4);
                delayTarget{delayNumber}(row,:)=target(1:35,4);
                delayFoil{delayNumber}(row,:)=foil(1:35,4);
                delayRT{delayNumber}{row}=reinfTarget;
                delayRF{delayNumber}{row}=reinfFoil;
                delayT{delayNumber}{row}=target;
                delayF{delayNumber}{row}=foil;
                delayCounter(delayNumber)=row+1;
            end
        end
    end

    numberOfMice=size(optomeanMat,1)-1;
    for delayNumber=1:5
        sessionKeys=delaySessionKey{delayNumber};
        expectedRows=0;
        for mouseIdx=1:numberOfMice
            for sessionId=1:size(delayBySession,1)
                expectedRows=expectedRows+nnz(delayBySession{sessionId,mouseIdx}==delayNumber);
            end
        end
        assert(size(sessionKeys,1)==expectedRows, ...
            'Missing session rows for delay %d: expected %d, found %d.', ...
            delayNumber,expectedRows,size(sessionKeys,1));
        assert(size(unique(sessionKeys,'rows'),1)==size(sessionKeys,1), ...
            'Duplicate mouse/session rows found for delay %d.',delayNumber);
    end
    allMiceDelay1RTargetIdx=delayRTargetIdx{1}; allMiceDelay1RFoilIdx=delayRFoilIdx{1};
    allMiceDelay1TargetIdx=delayTargetIdx{1}; allMiceDelay1FoilIdx=delayFoilIdx{1};
    allMiceDelay1RTarget=delayRTarget{1}; allMiceDelay1RFoil=delayRFoil{1};
    allMiceDelay1Target=delayTarget{1}; allMiceDelay1Foil=delayFoil{1};
    allMiceDelay1RT=delayRT{1}; allMiceDelay1RF=delayRF{1};
    allMiceDelay1T=delayT{1}; allMiceDelay1F=delayF{1};

    allMiceDelay2RTargetIdx=delayRTargetIdx{2}; allMiceDelay2RFoilIdx=delayRFoilIdx{2};
    allMiceDelay2TargetIdx=delayTargetIdx{2}; allMiceDelay2FoilIdx=delayFoilIdx{2};
    allMiceDelay2RTarget=delayRTarget{2}; allMiceDelay2RFoil=delayRFoil{2};
    allMiceDelay2Target=delayTarget{2}; allMiceDelay2Foil=delayFoil{2};
    allMiceDelay2RT=delayRT{2}; allMiceDelay2RF=delayRF{2};
    allMiceDelay2T=delayT{2}; allMiceDelay2F=delayF{2};

    allMiceDelay3RTargetIdx=delayRTargetIdx{3}; allMiceDelay3RFoilIdx=delayRFoilIdx{3};
    allMiceDelay3TargetIdx=delayTargetIdx{3}; allMiceDelay3FoilIdx=delayFoilIdx{3};
    allMiceDelay3RTarget=delayRTarget{3}; allMiceDelay3RFoil=delayRFoil{3};
    allMiceDelay3Target=delayTarget{3}; allMiceDelay3Foil=delayFoil{3};
    allMiceDelay3RT=delayRT{3}; allMiceDelay3RF=delayRF{3};
    allMiceDelay3T=delayT{3}; allMiceDelay3F=delayF{3};

    allMiceDelay4RTargetIdx=delayRTargetIdx{4}; allMiceDelay4RFoilIdx=delayRFoilIdx{4};
    allMiceDelay4TargetIdx=delayTargetIdx{4}; allMiceDelay4FoilIdx=delayFoilIdx{4};
    allMiceDelay4RTarget=delayRTarget{4}; allMiceDelay4RFoil=delayRFoil{4};
    allMiceDelay4Target=delayTarget{4}; allMiceDelay4Foil=delayFoil{4};
    allMiceDelay4RT=delayRT{4}; allMiceDelay4RF=delayRF{4};
    allMiceDelay4T=delayT{4}; allMiceDelay4F=delayF{4};

    allMiceDelay5RTargetIdx=delayRTargetIdx{5}; allMiceDelay5RFoilIdx=delayRFoilIdx{5};
    allMiceDelay5TargetIdx=delayTargetIdx{5}; allMiceDelay5FoilIdx=delayFoilIdx{5};
    allMiceDelay5RTarget=delayRTarget{5}; allMiceDelay5RFoil=delayRFoil{5};
    allMiceDelay5Target=delayTarget{5}; allMiceDelay5Foil=delayFoil{5};
    allMiceDelay5RT=delayRT{5}; allMiceDelay5RF=delayRF{5};
    allMiceDelay5T=delayT{5}; allMiceDelay5F=delayF{5};

SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% all animals, by individual, scatterplt %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %first reorganinze the lick mat file 
    count=0;
    for qq=1:length(allLicksTest)
        for aa=1:length(allLicksTest{1,qq})
            count=count+1;
            if aa==1
                allLicks{count,1}=allLicksTest{1,qq}{aa,1};
                allLicks{count,2}=allLicksTest{1,qq}{aa,2};
            else
                temp=allLicksTest{1,qq}{aa,1};
                allLicks{count,1}=temp(2:end,:);
                temp=allLicksTest{1,qq}{aa,2};
                allLicks{count,2}=temp;
            end
        end
    end
    allLicksByConditionVert=vertcat(allLicks{:,1}); % this is by each session (10) and animal
    allLicksVert=vertcat(allLicks{:,2});
    count=0;
    allDays={'sk198','sk203','sk204';[0,5,1],[0,2,1],[0,2,1];...
        [3,4,0],[3,4,0],[3,4,0];...
        [5,1,0],[5,1,0],[5,2,0];[5,0,2],[5,2,0],[5,1,0];...
        [3,4,0],[3,4,0],[3,0,4]};
    ctxtOrder=[1,5,6];
    % allDays is organized by session and opto condition (ctxt 0, ctxt 5,ctxt 6). 
    % '0,2,1' means for that session, there were 
    % no full trial conditions (context 1), 
    % delay 2 was tone (context 5)
    % delay 1 was choice (context 6)
    for nbsubj=1:3
        [tempLickMat,rD1TLickMat,rD1FLickMat,D1TLickMat,D1FLickMat,...
        rD2TLickMat, rD2FLickMat, D2TLickMat, D2FLickMat, rD3TLickMat, rD3FLickMat, ...
        D3TLickMat, D3FLickMat, rD4TLickMat, rD4FLickMat, D4TLickMat, D4FLickMat, ...
        rD5TLickMat, rD5FLickMat, D5TLickMat, D5FLickMat]=initDelayMat;

        days=allDays(2:6,nbsubj);
        subjDays=vertcat(days{:});
        count=count+1;
        % now get the relevant days
        [delay1Day,~]=find(subjDays==1); 
        [delay2Day,~]=find(subjDays==2); 
        [delay3Day,~]=find(subjDays==3); 
        [delay4Day,~]=find(subjDays==4); 
        [delay5Day,~]=find(subjDays==5); 
        if nbsubj==1
            sessRange=(1:length(days))+nbsubj;
            countday=4;
        else
            sessRange=(1:length(days))+nbsubj+countday;
            countday=countday+4;
        end
        STIM=5;CHOICE=7;FULL=3;
        animalConditionLicks=allLicksByConditionVert(sessRange(1):sessRange(5),:);
        animalAllLicks=allLicksVert((sessRange(1):sessRange(5))-1,:);
        for gg=1:size(subjDays,1) %make this flexible to iterate through all delay types
            if any(subjDays(gg,:)==1) % delay 1
                conditions=allDays{gg+1,nbsubj};
                optoCtxt=find(conditions==1);
                optoConditionbyCtxt=ctxtOrder(optoCtxt);
                licksRD1T=animalConditionLicks{gg,1};
                licksRD1F=animalConditionLicks{gg,2};
                if optoCtxt==2 % this is a stimulus condition
                    licksD1T=animalConditionLicks{gg,STIM};
                    licksD1F=animalConditionLicks{gg,STIM+1};
                    ctxtFlag=STIM;
                elseif optoCtxt==1
                    % this is a full condition
                    licksD1T=animalConditionLicks{gg,FULL};
                    licksD1F=animalConditionLicks{gg,FULL+1};
                    ctxtFlag=FULL;
                elseif optoCtxt==3 % choice condition
                    licksD1T=animalConditionLicks{gg,CHOICE};
                    licksD1F=animalConditionLicks{gg,CHOICE+1};
                    ctxtFlag=CHOICE;
                end

            elseif any(subjDays(gg,:)==2) % delay 2
                conditions=allDays{gg+1,nbsubj};
                optoCtxt=find(conditions==2);
                optoConditionbyCtxt=ctxtOrder(optoCtxt);
                licksRD2T=animalConditionLicks{gg,1};
                licksRD2F=animalConditionLicks{gg,2};
                if optoCtxt==2 % this is a stimulus condition
                    ctxtFlag=STIM;
                    licksD2T=animalConditionLicks{gg,STIM};
                    licksD2F=animalConditionLicks{gg,STIM+1};
                elseif optoCtxt==1
                    % this is a full condition
                    ctxtFlag=FULL;
                    licksD2T=animalConditionLicks{gg,FULL};
                    licksD2F=animalConditionLicks{gg,FULL+1};
                elseif optoCtxt==3 % choice condition
                    ctxtFlag=CHOICE;
                    licksD2T=animalConditionLicks{gg,CHOICE};
                    licksD2F=animalConditionLicks{gg,CHOICE+1};
                end
            elseif any(subjDays(gg,:)==3) % delay 3
                conditions=allDays{gg+1,nbsubj};
                optoCtxt=find(conditions==3);
                optoConditionbyCtxt=ctxtOrder(optoCtxt);
                licksRD3T=animalConditionLicks{gg,1};
                licksRD3F=animalConditionLicks{gg,2};
                if optoCtxt==2 % this is a stimulus condition
                    ctxtFlag=STIM;
                    licksD3T=animalConditionLicks{gg,STIM};
                    licksD3F=animalConditionLicks{gg,STIM+1};
                elseif optoCtxt==1
                    % this is a full condition
                    ctxtFlag=FULL;
                    licksD3T=animalConditionLicks{gg,FULL};
                    licksD3F=animalConditionLicks{gg,FULL+1};
                elseif optoCtxt==3 % choice condition
                    ctxtFlag=CHOICE;
                    licksD3T=animalConditionLicks{gg,CHOICE};
                    licksD3F=animalConditionLicks{gg,CHOICE+1};
                end
            elseif any(subjDays(gg,:)==4) % delay 4
                conditions=allDays{gg+1,nbsubj};
                optoCtxt=find(conditions==4);
                optoConditionbyCtxt=ctxtOrder(optoCtxt);
                licksRD4T=animalConditionLicks{gg,1};
                licksRD4F=animalConditionLicks{gg,2};
                if optoCtxt==2 % this is a stimulus condition
                    ctxtFlag=STIM;
                    licksD4T=animalConditionLicks{gg,STIM};
                    licksD4F=animalConditionLicks{gg,STIM+1};
                elseif optoCtxt==1
                    % this is a full condition
                    ctxtFlag=FULL;
                    licksD4T=animalConditionLicks{gg,FULL};
                    licksD4F=animalConditionLicks{gg,FULL+1};
                elseif optoCtxt==3 % choice condition
                    ctxtFlag=CHOICE;
                    licksD4T=animalConditionLicks{gg,CHOICE};
                    licksD4F=animalConditionLicks{gg,CHOICE+1};
                end
            elseif any(subjDays(gg,:)==5) % delay 5
                conditions=allDays{gg+1,nbsubj};
                optoCtxt=find(conditions==5);
                optoConditionbyCtxt=ctxtOrder(optoCtxt);
                licksRD5T=animalConditionLicks{gg,1};
                licksRD5F=animalConditionLicks{gg,2};
                if optoCtxt==2 % this is a stimulus condition
                    ctxtFlag=STIM;licksD5T=animalConditionLicks{gg,STIM};
                    licksD5F=animalConditionLicks{gg,STIM+1};
                elseif optoCtxt==1
                    % this is a full condition
                    ctxtFlag=FULL;
                    licksD5T=animalConditionLicks{gg,FULL};
                    licksD5F=animalConditionLicks{gg,FULL+1};
                elseif optoCtxt==3 % choice condition
                     ctxtFlag=CHOICE;
                    licksD5T=animalConditionLicks{gg,CHOICE};
                    licksD5F=animalConditionLicks{gg,CHOICE+1};
                end
            end
            
            % get consummatory licks
            sessLicks=animalAllLicks{gg};
            exampleSession=optomeanMat{nbsubj+1,17}; 
            exampleTrials=find(exampleSession(:,1)==gg);
            exampleTrials=exampleSession(exampleTrials,:);
            
%             nextIdx=1; 
%             for tt=1:size(exampleTrials)
%                 nextIdxTemp=find(sessLicks>exampleTrials(tt,6));
%                 try
%                     nextIdx(tt+1)=nextIdxTemp(1);
%                 catch
%                     disp(tt)
%                 end
%             end

            nTrials = size(exampleTrials,1);
            nextIdx = nan(1,nTrials);
            nextIdx(1) = 1;

            for tt = 1:nTrials
                idx = find(sessLicks > exampleTrials(tt,6), 1, 'first');

                if isempty(idx)
                    % End boundary: one position beyond the final lick
                    nextIdx(tt+1) = numel(sessLicks) + 1;

                    warning('No lick after trial %d; trial time=%g, last lick=%g', ...
                        tt, exampleTrials(tt,6), sessLicks(end));
                else
                    nextIdx(tt+1) = idx;
                end
            end
            
            for yu=1:length(nextIdx)-1
                tempLickMat(yu,1:length(sessLicks(nextIdx(yu):nextIdx(yu+1)-1)))=sessLicks(nextIdx(yu):nextIdx(yu+1)-1);
            end

            %now tempLickMat is the licks, for each trial, for the entire session
            % sort tempLickMat now by tone and context/condition
            if any(subjDays(gg,:)==1) % delay 1
                reinfD1TIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==1);
                reinfD1FIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==2);
                reinfD1TLicks=tempLickMat(reinfD1TIdx,:);
                reinfD1FLicks=tempLickMat(reinfD1FIdx,:);
                for ym=1:length(nextIdx)-1
                    rD1TLickMat(ym,1:length(licksRD1T(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD1T(nextIdx(ym):nextIdx(ym+1)-1);
                    rD1FLickMat(ym,1:length(licksRD1F(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD1F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                rD1TLickMatAll{gg}=rD1TLickMat(reinfD1TIdx,1:60); %there's a bug somewhere making a ton of 0's--where?
                rD1FLickMatAll{gg}=rD1FLickMat(reinfD1FIdx,1:60);

                D1TIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==1);
                D1FIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==2);
                D1TLicks=tempLickMat(D1TIdx,:);
                D1FLicks=tempLickMat(D1FIdx,:);
                for ym=1:length(nextIdx)-1
                    D1TLickMat(ym,1:length(licksD1T(nextIdx(ym):nextIdx(ym+1)-1)))=licksD1T(nextIdx(ym):nextIdx(ym+1)-1);
                    D1FLickMat(ym,1:length(licksD1F(nextIdx(ym):nextIdx(ym+1)-1)))=licksD1F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                D1TLickMatAll{gg}=D1TLickMat(D1TIdx,1:60);
                D1FLickMatAll{gg}=D1FLickMat(D1FIdx,1:60);

                
            elseif any(subjDays(gg,:)==2) % delay 2
                reinfD2TIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==1);
                reinfD2FIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==2);
                reinfD2TLicks=tempLickMat(reinfD2TIdx,:);
                reinfD2FLicks=tempLickMat(reinfD2FIdx,:);
                for ym=1:length(nextIdx)-1
                    rD2TLickMat(ym,1:length(licksRD2T(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD2T(nextIdx(ym):nextIdx(ym+1)-1);
                    rD2FLickMat(ym,1:length(licksRD2F(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD2F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                rD2TLickMatAll{gg}=rD2TLickMat(reinfD2TIdx,1:60);
                rD2FLickMatAll{gg}=rD2FLickMat(reinfD2FIdx,1:60);

                D2TIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==1);
                D2FIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==2);
                D2TLicks=tempLickMat(D2TIdx,:);
                D2FLicks=tempLickMat(D2FIdx,:);
                for ym=1:length(nextIdx)-1
                    D2TLickMat(ym,1:length(licksD2T(nextIdx(ym):nextIdx(ym+1)-1)))=licksD2T(nextIdx(ym):nextIdx(ym+1)-1);
                    D2FLickMat(ym,1:length(licksD2F(nextIdx(ym):nextIdx(ym+1)-1)))=licksD2F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                D2TLickMatAll{gg}=D2TLickMat(D2TIdx,1:60);
                D2FLickMatAll{gg}=D2FLickMat(D2FIdx,1:60);
            elseif any(subjDays(gg,:)==3) % delay 3
                reinfD3TIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==1);
                reinfD3FIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==2);
                reinfD3TLicks=tempLickMat(reinfD3TIdx,:);
                reinfD3FLicks=tempLickMat(reinfD3FIdx,:);
                for ym=1:length(nextIdx)-1
                    rD3TLickMat(ym,1:length(licksRD3T(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD3T(nextIdx(ym):nextIdx(ym+1)-1);
                    rD3FLickMat(ym,1:length(licksRD3F(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD3F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                rD3TLickMatAll{gg}=rD3TLickMat(reinfD3TIdx,1:60);
                rD3FLickMatAll{gg}=rD3FLickMat(reinfD3FIdx,1:60);

                D3TIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==1);
                D3FIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==2);
                D3TLicks=tempLickMat(D3TIdx,:);
                D3FLicks=tempLickMat(D3FIdx,:);
                for ym=1:length(nextIdx)-1
                    D3TLickMat(ym,1:length(licksD3T(nextIdx(ym):nextIdx(ym+1)-1)))=licksD3T(nextIdx(ym):nextIdx(ym+1)-1);
                    D3FLickMat(ym,1:length(licksD3F(nextIdx(ym):nextIdx(ym+1)-1)))=licksD3F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                D3TLickMatAll{gg}=D3TLickMat(D3TIdx,1:60);
                D3FLickMatAll{gg}=D3FLickMat(D3FIdx,1:60);
            elseif any(subjDays(gg,:)==4) % delay 4
                reinfD4TIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==1);
                reinfD4FIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==2);
                reinfD4TLicks=tempLickMat(reinfD4TIdx,:);
                reinfD4FLicks=tempLickMat(reinfD4FIdx,:);
                for ym=1:length(nextIdx)-1
                    rD4TLickMat(ym,1:length(licksRD4T(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD4T(nextIdx(ym):nextIdx(ym+1)-1);
                    rD4FLickMat(ym,1:length(licksRD4F(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD4F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                rD4TLickMat=rD4TLickMat(reinfD4TIdx,1:60);
                rD4FLickMat=rD4FLickMat(reinfD4FIdx,1:60);

                D4TIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==1);
                D4FIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==2);
                D4TLicks=tempLickMat(D4TIdx,:);
                D4FLicks=tempLickMat(D4FIdx,:);
                for ym=1:length(nextIdx)-1
                    D4TLickMat(ym,1:length(licksD4T(nextIdx(ym):nextIdx(ym+1)-1)))=licksD4T(nextIdx(ym):nextIdx(ym+1)-1);
                    D4FLickMat(ym,1:length(licksD4F(nextIdx(ym):nextIdx(ym+1)-1)))=licksD4F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                D4TLickMatAll(gg)=D4TLickMat(D4TIdx,1:60);
                D4FLickMatAll(gg)=D4FLickMat(D4FIdx,1:60);
            elseif any(subjDays(gg,:)==5) % delay 5
                reinfD5TIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==1);
                reinfD5FIdx=find(exampleTrials(1:300,CTXT)==2 & exampleTrials(1:300,3)==2);
                reinfD5TLicks=tempLickMat(reinfD5TIdx,:);
                reinfD5FLicks=tempLickMat(reinfD5FIdx,:);
                for ym=1:length(nextIdx)-1
                    rD5TLickMat(ym,1:length(licksRD5T(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD5T(nextIdx(ym):nextIdx(ym+1)-1);
                    rD5FLickMat(ym,1:length(licksRD5F(nextIdx(ym):nextIdx(ym+1)-1)))=licksRD5F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                rD5TLickMatAll{gg}=rD5TLickMat(reinfD5TIdx,1:60);
                rD5FLickMatAll{gg}=rD5FLickMat(reinfD5FIdx,1:60);

                D5TIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==1);
                D5FIdx=find(exampleTrials(1:300,CTXT)==optoConditionbyCtxt & exampleTrials(1:300,3)==2);
                D5TLicks=tempLickMat(D5TIdx,:);
                D5FLicks=tempLickMat(D5FIdx,:);
                for ym=1:length(nextIdx)-1
                    D5TLickMat(ym,1:length(licksD5T(nextIdx(ym):nextIdx(ym+1)-1)))=licksD5T(nextIdx(ym):nextIdx(ym+1)-1);
                    D5FLickMat(ym,1:length(licksD5F(nextIdx(ym):nextIdx(ym+1)-1)))=licksD5F(nextIdx(ym):nextIdx(ym+1)-1);
                end
                D5TLickMatAll(gg)=D5TLickMat(D5TIdx,1:60);
                D5FLickMatAll(gg)=D5FLickMat(D5FIdx,1:60);
            end
            
        end

        % this part does not need to be in the loop
        % this is fixed 7/3/26
        lickLatRD1T=allMiceDelay1RT{1,count}(:,LICKL);
        lickLatRD1F=allMiceDelay1RF{1,count}(:,LICKL);
        lickLatD1T=allMiceDelay1T{1,count}(:,LICKL);
        lickLatD1F=allMiceDelay1F{1,count}(:,LICKL);
        lickLatRD1TNoNan=lickLatRD1T(~isnan(lickLatRD1T));
        lickLatRD1FNoNan=lickLatRD1F(~isnan(lickLatRD1F));
        lickLatD1TNoNan=lickLatD1T(~isnan(lickLatD1T));
        lickLatD1FNoNan=lickLatD1F(~isnan(lickLatD1F));
        
        lickLatRD2T=allMiceDelay2RT{1,count}(:,LICKL);
        lickLatRD2F=allMiceDelay2RF{1,count}(:,LICKL);
        lickLatD2T=allMiceDelay2T{1,count}(:,LICKL);
        lickLatD2F=allMiceDelay2F{1,count}(:,LICKL);
        lickLatRD2TNoNan=lickLatRD2T(~isnan(lickLatRD2T));
        lickLatRD2FNoNan=lickLatRD2F(~isnan(lickLatRD2F));
        lickLatD2TNoNan=lickLatD2T(~isnan(lickLatD2T));
        lickLatD2FNoNan=lickLatD2F(~isnan(lickLatD2F));

        lickLatRD3T=allMiceDelay3RT{1,count}(:,LICKL);
        lickLatRD3F=allMiceDelay3RF{1,count}(:,LICKL);
        lickLatD3T=allMiceDelay3T{1,count}(:,LICKL);
        lickLatD3F=allMiceDelay3F{1,count}(:,LICKL);
        lickLatRD3TNoNan=lickLatRD3T(~isnan(lickLatRD3T));
        lickLatRD3FNoNan=lickLatRD3F(~isnan(lickLatRD3F));
        lickLatD3TNoNan=lickLatD3T(~isnan(lickLatD3T));
        lickLatD3FNoNan=lickLatD3F(~isnan(lickLatD3F));

        lickLatRD4T=allMiceDelay4RT{1,count}(:,LICKL);
        lickLatRD4F=allMiceDelay4RF{1,count}(:,LICKL);
        lickLatD4T=allMiceDelay4T{1,count}(:,LICKL);
        lickLatD4F=allMiceDelay4F{1,count}(:,LICKL);
        lickLatRD4TNoNan=lickLatRD4T(~isnan(lickLatRD4T));
        lickLatRD4FNoNan=lickLatRD4F(~isnan(lickLatRD4F));
        lickLatD4TNoNan=lickLatD4T(~isnan(lickLatD4T));
        lickLatD4FNoNan=lickLatD4F(~isnan(lickLatD4F));

        lickLatRD5T=allMiceDelay5RT{1,count}(:,LICKL);
        lickLatRD5F=allMiceDelay5RF{1,count}(:,LICKL);
        lickLatD5T=allMiceDelay5T{1,count}(:,LICKL);
        lickLatD5F=allMiceDelay5F{1,count}(:,LICKL);
        lickLatRD5TNoNan=lickLatRD5T(~isnan(lickLatRD5T));
        lickLatRD5FNoNan=lickLatRD5F(~isnan(lickLatRD5F));
        lickLatD5TNoNan=lickLatD5T(~isnan(lickLatD5T));
        lickLatD5FNoNan=lickLatD5F(~isnan(lickLatD5F));

        % now remove the probe block indices
        %'Condition','Trial indices','No probe indicies','Miss idx no probe'
        idxAllMice={'Reinf Delay1 T',allMiceDelay1RTargetIdx;'Reinf Delay1 F',allMiceDelay1RFoilIdx;...
            'Delay1 T',allMiceDelay1TargetIdx;'Delay1 F',allMiceDelay1FoilIdx;...
            'Reinf Delay2 T',allMiceDelay2RTargetIdx;'Reinf Delay2 F',allMiceDelay2RFoilIdx;...
            'Delay2 T',allMiceDelay2TargetIdx;'Delay2 F',allMiceDelay2FoilIdx;...
            'Reinf Delay3 T',allMiceDelay3RTargetIdx;'Reinf Delay3 F',allMiceDelay3RFoilIdx;...
            'Delay3 T',allMiceDelay3TargetIdx;'Delay3 F',allMiceDelay3FoilIdx;...
            'Reinf Delay4 T',allMiceDelay4RTargetIdx;'Reinf Delay4 F',allMiceDelay4RFoilIdx;...
            'Delay4 T',allMiceDelay4TargetIdx;'Delay4 F',allMiceDelay4FoilIdx;...
            'Reinf Delay5 T',allMiceDelay5RTargetIdx;'Reinf Delay5 F',allMiceDelay5RFoilIdx;...
            'Delay5 T',allMiceDelay5TargetIdx;'Delay5 F',allMiceDelay5FoilIdx};

        % column 4 records [mouse index, session index] for every matrix row
        for delayNumber=1:5
            rowOffset=(delayNumber-1)*4;
            for conditionRow=1:4
                idxAllMice{rowOffset+conditionRow,4}=delaySessionKey{delayNumber};
            end
        end
        
            for yi=1:length(idxAllMice)
                for oo=1:size(idxAllMice{yi,2},1)
                    probeIdx=find(idxAllMice{yi,2}(oo,:)>140);
                    idxAllMice{yi,3}(oo,:)=[idxAllMice{yi,2}(oo,1:probeIdx(1)-1) ...
                        (idxAllMice{yi,2}(oo,probeIdx(1):length(idxAllMice{yi,2}))-20)];
                end
            end


            %add open circles aligned to 0 for the misses
            missRD1Idx=find(allMiceDelay1RTarget(count,:)==0);
            missRD2Idx=find(allMiceDelay2RTarget(count,:)==0);
            missRD3Idx=find(allMiceDelay3RTarget(count,:)==0);
            missRD4Idx=find(allMiceDelay4RTarget(count,:)==0);
            missRD5Idx=find(allMiceDelay5RTarget(count,:)==0);
            
            missD1Idx=find(allMiceDelay1Target(count,:)==0);
            missD2Idx=find(allMiceDelay2Target(count,:)==0);
            missD3Idx=find(allMiceDelay3Target(count,:)==0);
            missD4dx=find(allMiceDelay4Target(count,:)==0);
            missD5Idx=find(allMiceDelay5Target(count,:)==0);
            
            crRD1Idx=find(allMiceDelay1RFoil(count,:)==0);
            crRD2Idx=find(allMiceDelay2RFoil(count,:)==0);
            crRD3Idx=find(allMiceDelay3RFoil(count,:)==0);
            crRD4Idx=find(allMiceDelay4RFoil(count,:)==0);
            crRD5Idx=find(allMiceDelay5RFoil(count,:)==0);
            
            crD1Idx=find(allMiceDelay1Foil(count,:)==0);
            crD2Idx=find(allMiceDelay2Foil(count,:)==0);
            crD3Idx=find(allMiceDelay3Foil(count,:)==0);
            crD4Idx=find(allMiceDelay4Foil(count,:)==0);
            crD5Idx=find(allMiceDelay5Foil(count,:)==0);
            % need to fix this to use the new indexing without probe
            animalCell=allDays{1,nbsubj};
            sz=10;licksColor=[0.9 0.9 0.9];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %stopped here 7/7/26
            scatterPSTH=1;missColor=[1.0, 0.27, 0.0];
            if scatterPSTH==1
                hitFullFigD1=figure;
                delayLabel=find(cellfun(@(x) contains(x, 'Delay1'), idxAllMice(:,1)));
                rCount=idxAllMice{delayLabel(1),4}(count,2);
                subplot(2,3,1);scatter(rD1TLickMatAll{count},idxAllMice{delayLabel(1),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatRD1T,idxAllMice{delayLabel(1),3}(count,:)',sz,reinfcolor,'filled');
                title('Light off, Hit');xlim([-0.5 3]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missRD1Idx))*3,idxAllMice{delayLabel(1),3}(count,missRD1Idx)',sz,missColor);
                subplot(2,3,2);scatter(D1TLickMatAll{count},idxAllMice{delayLabel(3),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatD1T,idxAllMice{delayLabel(3),3}(count,:)',sz,optocolor,'filled'); title('Delay 1, Hit');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missD1Idx))*3,idxAllMice{delayLabel(3),3}(count,missD1Idx)',sz,missColor);ylabel('Trial');
                subplot(2,3,3);shadedErrorBar(1:length(lickLatRD1TNoNan),lickLatRD1TNoNan,std(lickLatRD1TNoNan)); 
                hold on;shadedErrorBar(1:length(lickLatD1TNoNan),lickLatD1TNoNan(1:length(lickLatD1TNoNan)),std(lickLatD1TNoNan),'b');ylabel('First lick latency (s)');
                xlim([1 20]);box off;title('Delay 1, Hit');
                subplot(2,3,4);scatter(rD1FLickMatAll{count},idxAllMice{delayLabel(2),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatRD1F,idxAllMice{delayLabel(2),3}(count,:)',sz,reinfcolor,'filled');title('Light off, FA');xlim([0 2]);
                ylabel('Trial');xlabel('Time (s)');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(crRD1Idx))*3,idxAllMice{2,3}(count,crRD1Idx)',sz,reinfcolor);ylabel('Trial');
                subplot(2,3,5);scatter(D1FLickMatAll{count},idxAllMice{delayLabel(4),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatD1F,idxAllMice{delayLabel(4),3}(count,:)',sz,optocolor,'filled'); title('Delay 1, FA');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missD1Idx))*3,idxAllMice{delayLabel(4),3}(count,missD1Idx)',sz,optocolor);ylabel('Trial');
                subplot(2,3,6);shadedErrorBar(1:length(lickLatRD1FNoNan),lickLatRD1FNoNan,std(lickLatRD1FNoNan)); 
                hold on;shadedErrorBar(1:length(lickLatD1FNoNan),lickLatD1FNoNan(1:length(lickLatD1FNoNan)),std(lickLatD1FNoNan),'b');ylabel('First lick latency (s)');
                xlim([1 20]);box off;title('Delay 1, FA');
                hitFullFigD1.Position(3:4)=[550 350];
                saveas(gcf,[animalCell 'D' num2str(gg) '_T_MGB_Delay1']);
                saveas(gcf,[animalCell 'D' num2str(gg) '_T_MGB_Delay1.png']);    

                %%%%%%%%%%%% Delay 2
                FigD2=figure;
                delayLabel=find(cellfun(@(x) contains(x, 'Delay2'), idxAllMice(:,1)));
                rCount=idxAllMice{delayLabel(1),4}(count,2);
                subplot(2,3,1);scatter(rD2TLickMatAll{rCount},idxAllMice{delayLabel(1),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatRD2T,idxAllMice{delayLabel(1),3}(count,:)',sz,reinfcolor,'filled');
                title('Light off, Hit');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missRD2Idx))*3,idxAllMice{delayLabel(1),3}(count,missRD2Idx)',sz,missColor);
                subplot(2,3,2);scatter(D2TLickMatAll{rCount},idxAllMice{delayLabel(3),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatD2T,idxAllMice{delayLabel(3),3}(count,:)',sz,optocolor,'filled'); title('Delay 2, Hit');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missD2Idx))*3,idxAllMice{delayLabel(3),3}(count,missD2Idx)',sz,missColor);ylabel('Trial');
                subplot(2,3,3);shadedErrorBar(1:length(lickLatRD2TNoNan),lickLatRD2TNoNan,std(lickLatRD2TNoNan)); 
                hold on;shadedErrorBar(1:length(lickLatD2TNoNan),lickLatD2TNoNan(1:length(lickLatD2TNoNan)),std(lickLatD2TNoNan),'b');ylabel('First lick latency (s)');
                xlim([1 20]);box off;title('Delay 2, Hit');
                subplot(2,3,4);scatter(rD2FLickMatAll{rCount},idxAllMice{delayLabel(2),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatRD2F,idxAllMice{delayLabel(2),3}(count,:)',sz,reinfcolor,'filled');title('Light off, FA');xlim([0 2]);
                ylabel('Trial');xlabel('Time (s)');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(crRD2Idx))*3,idxAllMice{delayLabel(2),3}(count,crRD2Idx)',sz,optocolor);ylabel('Trial');
                subplot(2,3,5);scatter(D2FLickMatAll{rCount},idxAllMice{delayLabel(4),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatD2F,idxAllMice{delayLabel(4),3}(count,:)',sz,optocolor,'filled'); title('Delay 2, FA');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missD2Idx))*3,idxAllMice{delayLabel(4),3}(count,missD2Idx)',sz,optocolor);ylabel('Trial');
                subplot(2,3,6);shadedErrorBar(1:length(lickLatRD2FNoNan),lickLatRD2FNoNan,std(lickLatRD2FNoNan)); 
                hold on;shadedErrorBar(1:length(lickLatD2FNoNan),lickLatD2FNoNan(1:length(lickLatD2FNoNan)),std(lickLatD2FNoNan),'b');ylabel('First lick latency (s)');
                xlim([1 20]);box off;title('Delay 2, FA');
                FigD2.Position(3:4)=[550 350];
                saveas(gcf,[animalCell 'D' num2str(gg) '_Delay2']);
                saveas(gcf,[animalCell 'D' num2str(gg) '_Delay2.png']);

                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % left off here 7/17/26
                %%%%%%%%%%%% Delay 3
                FigD3=figure;
                delayLabel=find(cellfun(@(x) contains(x, 'Delay3'), idxAllMice(:,1)));
                rCount=idxAllMice{1,4}(delayLabel(1),2);
                subplot(2,3,1);scatter(rD2TLickMatAll{rCount},idxAllMice{delayLabel(1),3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatRD2T,idxAllMice{delayLabel(1),3}(count,:)',sz,reinfcolor,'filled');
                title('Light off, Hit');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missRD2Idx))*3,idxAllMice{5,3}(count,missRD2Idx)',sz,missColor);

                subplot(2,3,2);scatter(D2TLickMatAll{rCount},idxAllMice{7,3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatD2T,idxAllMice{7,3}(count,:)',sz,optocolor,'filled'); title('Delay 2, Hit');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missD2Idx))*3,idxAllMice{7,3}(count,missD2Idx)',sz,missColor);ylabel('Trial');

                subplot(2,3,3);shadedErrorBar(1:length(lickLatRD2TNoNan),lickLatRD2TNoNan,std(lickLatRD2TNoNan)); 
                hold on;shadedErrorBar(1:length(lickLatD2TNoNan),lickLatD2TNoNan(1:length(lickLatD2TNoNan)),std(lickLatD2TNoNan),'b');ylabel('First lick latency (s)');
                xlim([1 20]);box off;title('Delay 2, Hit');

                subplot(2,3,4);scatter(rD2FLickMatAll{rCount},idxAllMice{6,3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatRD2F,idxAllMice{6,3}(count,:)',sz,reinfcolor,'filled');title('Light off, FA');xlim([0 2]);
                ylabel('Trial');xlabel('Time (s)');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(crRD2Idx))*3,idxAllMice{6,3}(count,crRD2Idx)',sz,optocolor);ylabel('Trial');

                subplot(2,3,5);scatter(D2FLickMatAll{rCount},idxAllMice{8,3}(count,:)',sz,licksColor,'filled');
                hold on;scatter(lickLatD2F,idxAllMice{8,3}(count,:)',sz,optocolor,'filled'); title('Delay 2, FA');xlim([0 2]);
                ylabel('Trial');xlim([-0.5 3]);ylim([0 280]);
                scatter(ones(1,length(missD2Idx))*3,idxAllMice{8,3}(count,missD2Idx)',sz,optocolor);ylabel('Trial');

                subplot(2,3,6);shadedErrorBar(1:length(lickLatRD2FNoNan),lickLatRD2FNoNan,std(lickLatRD2FNoNan)); 
                hold on;shadedErrorBar(1:length(lickLatD2FNoNan),lickLatD2FNoNan(1:length(lickLatD2FNoNan)),std(lickLatD2FNoNan),'b');ylabel('First lick latency (s)');
                xlim([1 20]);box off;title('Delay 2, FA');

                FigD2.Position(3:4)=[550 350];
                saveas(gcf,[animalCell 'D' num2str(gg) '_Delay2']);
                saveas(gcf,[animalCell 'D' num2str(gg) '_Delay2.png']);
            end
    end 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% all animals, lick latency scatterplt %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nbins=50;bins = linspace(-1,4,nbins);
    for pp=1:12
        lickLatRFT(1:70,pp)=allMiceRFT{1,pp}(:,LICKL);
        lickLatRFF(1:70,pp)=allMiceRFF{1,pp}(:,LICKL);
        lickLatFullT(1:35,pp)=allMiceFullT{1,pp}(:,LICKL);
        lickLatFullF(1:35,pp)=allMiceFullF{1,pp}(:,LICKL);
        lickLatCFullT(1:35,pp)=allMiceCFT{1,pp}(:,LICKL);
        lickLatCFullF(1:35,pp)=allMiceCFF{1,pp}(:,LICKL);
        
        lickLatRTT(1:70,pp)=allMiceRTT{1,pp}(:,LICKL);
        lickLatRTF(1:70,pp)=allMiceRTF{1,pp}(:,LICKL);
        lickLatToneT(1:35,pp)=allMiceToneT{1,pp}(:,LICKL);
        lickLatToneF(1:35,pp)=allMiceToneF{1,pp}(:,LICKL);
        lickLatCT(1:35,pp)=allMiceCTT{1,pp}(:,LICKL);
        lickLatCF(1:35,pp)=allMiceCTF{1,pp}(:,LICKL);
    end
    
    RFTHist=hist(lickLatRFT(:),bins)/(length(lickLatRFT)*6);
    RFFHist=hist(lickLatRFF(:),bins)/(length(lickLatRFF)*6);
    FullTHist=hist(lickLatFullT(:),bins)/(length(lickLatFullT)*6);
    FullFHist=hist(lickLatFullF(:),bins)/(length(lickLatFullF)*6);
    CFTHist=hist(lickLatCFullT(:),bins)/(length(lickLatCFullT)*6);
    CFFHist=hist(lickLatCFullF(:),bins)/(length(lickLatCFullF)*6);
       
    RTTHist=hist(lickLatRTT(:),bins)/(length(lickLatRTT)*6);
    RTFHist=hist(lickLatRTF(:),bins)/(length(lickLatRTF)*6);
    ToneTHist=hist(lickLatToneT(:),bins)/(length(lickLatToneT)*6);
    ToneFHist=hist(lickLatToneF(:),bins)/(length(lickLatToneF)*6);
    CTHist=hist(lickLatCT(:),bins)/(length(lickLatCT)*6);
    CFHist=hist(lickLatCF(:),bins)/(length(lickLatCF)*6);
    
  
    %add bar plot aligned to 4s for misses
    missRD1Idx=find(allMiceRFTarget==0);
    missFIdx=find(allMiceFullTarget==0);
    missFCIdx=find(allMiceCFTarget==0);
    missRTIdx=find(allMiceRTTarget==0);
    missTIdx=find(allMiceToneTarget==0);
    missTCIdx=find(allMiceCTTarget==0);
    crRD1Idx=find(allMiceRFFoil==0);
    crFIdx=find(allMiceFullFoil==0);
    crFCIdx=find(allMiceCFFoil==0);
    crRTIdx=find(allMiceRTFoil==0);
    crTIdx=find(allMiceToneFoil==0);
    crTCIdx=find(allMiceCTFoil==0);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% all animals, lick latency histograms & scatter plots %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    noLickColor=[0.8 0.3 0.3];
    lightOptoColor=[0.5,0.8,1];
    lighterOptoColor=[0.6,0.9,1];
    
    lickLatScatterHitFull=figure;
    subplot(2,3,1);title('Light off, Hit'); hold on;xlim([-0.5 4]);
    ylabel('Trial Number');xlabel('Time of first lick(s)');
    scatter(lickLatRFT,1:length(lickLatRFT),sz,reinfcolor);
    subplot(2,3,2); title('Full, Hit'); hold on;xlim([-0.5 4]);
    xlabel('Time of first lick(s)');ylim([0 35]);
    scatter(lickLatFullT,1:length(lickLatFullT),sz,optocolor);
    subplot(2,3,3); title('Choice, Hit'); hold on;xlabel('Time of first lick(s)');
    scatter(lickLatCFullT,1:length(lickLatCFullT),sz,optocolor);xlim([-0.5 4]);ylim([0 35]);
    subplot(2,3,4);
    plot(bins,RFTHist,'Color',reinfcolor,'LineWidth',2);xlim([-0.5 4]);
    ylabel('p(licks)');xlabel('Time of first lick(s)');ylim([0 0.6]);
    yyaxis right;ylabel('% of misses');
    hold on; rftMiss=bar(4,length(missRD1Idx)/(size(allMiceRFTarget,1)*size(allMiceRFTarget,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off;ylim([0 1]);
    subplot(2,3,5);hold off;
    plot(bins,FullTHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.6]);xlim([-0.5 4]);
    yyaxis right;ylabel('% of misses');
    hold on; rftMiss=bar(4,length(missFIdx)/(size(allMiceFullTarget,1)*size(allMiceFullTarget,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];ylim([0 1]);box off
    subplot(2,3,6);plot(bins,CFTHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.6]);xlim([-0.5 4]);
    yyaxis right;ylim([0 1]);
    hold on; rftMiss=bar(4,length(missFCIdx)/(size(allMiceCFTarget,1)*size(allMiceCFTarget,2)));box off
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];ylim([0 1]);
    ylabel('% of misses');
    lickLatScatterHitFull.Position(3:4)=[700 275];
    saveas(gcf,['allmice_T_MGB_HitFull_LickLatHist']);
    saveas(gcf,['allmice_T_MGB_HitFull_LickLatHist.png']);    
     saveas(gcf,['allmice_T_MGB_HitFull_LickLatHist.pdf']);  
     
    lickLatScatterHitTone=figure;
    subplot(2,3,1);title('Light off, Hit'); hold on;xlim([-0.5 4]);
    ylabel('Trial Number');xlabel('Time of first lick(s)');
    scatter(lickLatRTT,1:length(lickLatRTT),sz,reinfcolor);
    subplot(2,3,2); title('Stimulus, Hit'); hold on;xlim([-0.5 4]);
    xlabel('Time of first lick(s)');ylim([0 35]);
    scatter(lickLatToneT,1:length(lickLatToneT),sz,optocolor);
    subplot(2,3,3); title('Choice, Hit'); hold on;xlabel('Time of first lick(s)');
    scatter(lickLatCT,1:length(lickLatCT),sz,optocolor);xlim([-0.5 4]);ylim([0 35]);
    subplot(2,3,4);plot(bins,RTTHist,'Color',reinfcolor,'LineWidth',2);
    xlim([-0.5 4]);xlabel('Time of first lick(s)');ylabel('p(licks)');
    yyaxis right
    hold on; rftMiss=bar(4,length(missRTIdx)/(size(allMiceRTTarget,1)*size(allMiceRTTarget,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);ylabel('% of misses');
    subplot(2,3,5);plot(bins,ToneTHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.5]);xlim([-0.5 4]);
    ylim([0 0.6]);
    yyaxis right
    hold on; rftMiss=bar(4,length(missTIdx)/(size(allMiceToneTarget,1)*size(allMiceToneTarget,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);ylabel('% of misses');
    subplot(2,3,6);plot(bins,CTHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.5]);xlim([-0.5 4]);
    ylim([0 0.6]);
    yyaxis right
    hold on; rftMiss=bar(4,length(missTCIdx)/(size(allMiceCTTarget,1)*size(allMiceCTTarget,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);ylabel('% of misses');
    lickLatScatterHitTone.Position(3:4)=[700 275];
    saveas(gcf,['allmice_T_MGB_HitStim_LickLatHist']);
    saveas(gcf,['allmice_T_MGB_HitStim_LickLatHist.png']);  
    saveas(gcf,['allmice_T_MGB_HitStim_LickLatHist.pdf']);  
    
    lickLatScatterFAFull=figure;
    subplot(2,3,1);title('Light off, FA'); hold on;xlim([-0.5 4]);
    ylabel('Trial Number');xlabel('Time of first lick(s)');
    scatter(lickLatRFF,1:length(lickLatRFF),sz,reinfcolor);
    subplot(2,3,2); title('Full, FA'); hold on;xlim([-0.5 4]);
    xlabel('Time of first lick(s)');ylim([0 35]);
    scatter(lickLatFullF,1:length(lickLatFullF),sz,optocolor);ylim([0 35]);
    subplot(2,3,3); title('Choice, FA'); hold on;xlabel('Time of first lick(s)');
    scatter(lickLatCFullF,1:length(lickLatCFullF),sz,optocolor);xlim([-0.5 4]);ylim([0 35]);
    subplot(2,3,4);plot(bins,RFFHist,'Color',reinfcolor,'LineWidth',2);ylim([0 0.6]);
    ylabel('% of misses');
    xlim([-0.5 4]);xlabel('Time of first lick(s)');
    yyaxis right;ylabel('p(licks)');
    hold on; rftMiss=bar(4,length(crRD1Idx)/(size(allMiceRFFoil,1)*size(allMiceRFFoil,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);
    subplot(2,3,5);plot(bins,FullFHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.5]);xlim([-0.5 4]);
    ylim([0 0.6]);
    yyaxis right;
    hold on; rftMiss=bar(4,length(crFIdx)/(size(allMiceFullFoil,1)*size(allMiceFullFoil,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);
    ylabel('% of correct rejects');
    subplot(2,3,6);plot(bins,CFFHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.5]);xlim([-0.5 4]);
    ylim([0 0.6]);
    yyaxis right
    hold on; rftMiss=bar(4,length(crFCIdx)/(size(allMiceCFFoil,1)*size(allMiceCFFoil,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);
    ylabel('% of correct rejects');
    lickLatScatterFAFull.Position(3:4)=[700 275];
    saveas(gcf,['allmice_T_MGB_FAFull_LickLatHist']);
    saveas(gcf,['allmice_T_MGB_FAFull_LickLatHist.png']);    

    lickLatScatterFATone=figure;
    subplot(2,3,1);title('Light off, FA'); hold on;
    ylabel('Trial Number');xlabel('Time of first lick(s)');xlim([-0.5 4]);
    scatter(lickLatRTF,1:length(lickLatRTF),sz,reinfcolor);
    subplot(2,3,2); title('Stimulus, FA'); hold on;xlim([-0.5 4]);
    xlabel('Time of first lick(s)');ylim([0 35]);
    scatter(lickLatToneF,1:length(lickLatToneF),sz,optocolor);
    subplot(2,3,3); title('Choice, FA'); hold on;xlabel('Time of first lick(s)');
    scatter(lickLatCF,1:length(lickLatCF),sz,optocolor);xlim([-0.5 4]);ylim([0 35]);
    subplot(2,3,4);plot(bins,RTFHist,'Color',reinfcolor,'LineWidth',2);ylim([0 0.6]);ylabel('p(licks)');
    xlim([-0.5 4]);
    yyaxis right;ylim([0 1]);ylabel('% of correct rejects');
    hold on; rftMiss=bar(4,length(crRTIdx)/(size(allMiceRTFoil,1)*size(allMiceRTFoil,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);xlabel('Time of first lick(s)');
    subplot(2,3,5);plot(bins,ToneFHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.6]);xlim([-0.5 4]);
    yyaxis right;ylim([0 1]);ylabel('% of correct rejects');
    hold on; rftMiss=bar(4,length(crTIdx)/(size(allMiceToneFoil,1)*size(allMiceToneFoil,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);xlabel('Time of first lick(s)');
    
    subplot(2,3,6);plot(bins,CFHist,'Color',optocolor,'LineWidth',2);xlabel('Time of first lick(s)');ylim([0 0.6]);xlim([-0.5 4]);
    yyaxis right;ylim([0 1])
    hold on; rftMiss=bar(4,length(crTCIdx)/(size(allMiceCTFoil,1)*size(allMiceCTFoil,2)));
    rftMiss(1).FaceColor='flat'; rftMiss(1).CData=[noLickColor];box off
	ylim([0 1]);xlabel('Time of first lick(s)');ylabel('% of correct rejects');
    lickLatScatterFATone.Position(3:4)=[700 275];
    saveas(gcf,['allmice_T_MGB_FAStim_LickLatHist']);
    saveas(gcf,['allmice_T_MGB_FAStim_LickLatHist.png']);    
    saveas(gcf,['allmice_T_MGB_FAStim_LickLatHist.pdf']);  

    violinFig=figure;
    subplot(1,2,1);
    padlickLatFullT=nan(size(lickLatRFT,1),size(lickLatRFT,2));
    padlickLatFullT(1:size(lickLatFullT,1),1:size(lickLatFullT,2))=lickLatFullT;
    padlickLatToneT=nan(size(lickLatRFT,1),size(lickLatRFT,2));
    padlickLatToneT(1:size(lickLatToneT,1),1:size(lickLatToneT,2))=lickLatToneT;
    padlickLatCAllT=nan(size(lickLatRFT,1),size(lickLatRFT,2));
    padlickLatCAllT=[lickLatCFullT;lickLatCT];
    violins=violinplot([lickLatRFT(:),padlickLatFullT(:),padlickLatToneT(:),padlickLatCAllT(:)]);hold on;
    set(violins(1).ViolinPlot(:),'FaceColor',reinfcolor,'EdgeColor','none');
    set(violins(1).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',reinfcolor-0.2);
    violins(1).ViolinColor{:}=reinfcolor;
    darkOptoColor=optocolor-0.3;
    set(violins(2).ViolinPlot(:),'FaceColor',darkOptoColor,'EdgeColor','none');
    set(violins(2).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',optocolor-0.1);
    violins(2).ViolinColor{:}=optocolor-0.1;
    set(violins(3).ViolinPlot(:),'FaceColor',lightOptoColor,'EdgeColor','none');
    set(violins(3).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',lightOptoColor-0.2);
    violins(3).ViolinColor{:}=lightOptoColor;
    set(violins(4).ViolinPlot(:),'FaceColor',lighterOptoColor,'EdgeColor','none');
    set(violins(4).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',lighterOptoColor-0.1);
    violins(4).ViolinColor{:}=lighterOptoColor;xlim([0.5 4.5]);
    box off;
    ylabel('Time of first S+ lick (s)');xticklabels({'Light off','Full','Stimulus','Choice'});
    violinFig.Position(3:4)=[250 500];
    
    
    subplot(2,1,2);
    padlickLatFullF=nan(size(lickLatRFF,1),size(lickLatRFF,2));
    padlickLatFullF(1:size(lickLatFullF,1),1:size(lickLatFullF,2))=lickLatFullF;
    padlickLatToneF=nan(size(lickLatRFF,1),size(lickLatRFF,2));
    padlickLatToneF(1:size(lickLatToneF,1),1:size(lickLatToneF,2))=lickLatToneF;
    padlickLatCAllF=nan(size(lickLatRFF,1),size(lickLatRFF,2));
    padlickLatCAllF=[lickLatCFullF;lickLatCF];
    violins=violinplot([lickLatRFF(:),padlickLatFullF(:),padlickLatToneF(:),padlickLatCAllF(:)]);hold on;
    set(violins(1).ViolinPlot(:),'FaceColor',reinfcolor,'EdgeColor','none');
    set(violins(1).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',reinfcolor-0.2);
    violins(1).ViolinColor{:}=reinfcolor;
    darkOptoColor=optocolor-0.3;
    set(violins(2).ViolinPlot(:),'FaceColor',darkOptoColor,'EdgeColor','none');
    set(violins(2).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',optocolor-0.1);
    violins(2).ViolinColor{:}=optocolor-0.1;
    set(violins(3).ViolinPlot(:),'FaceColor',lightOptoColor,'EdgeColor','none');
    set(violins(3).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',lightOptoColor-0.2);
    violins(3).ViolinColor{:}=lightOptoColor;
    set(violins(4).ViolinPlot(:),'FaceColor',lighterOptoColor,'EdgeColor','none');
    set(violins(4).BoxPlot(:),'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'LineWidth',3,...
        'MarkerFaceColor',lighterOptoColor-0.1);
    violins(4).ViolinColor{:}=lighterOptoColor;xlim([0.5 4.5]);
    box off;
    ylabel('Time of first S+ lick (s)');xticklabels({'Light off','Full','Stimulus','Choice'});
    saveas(gcf,['allmice_T_MGB_Violin_LickLatSum']);
    saveas(gcf,['allmice_T_MGB_Violin_LickLatSum.png']);
        
    summaryLickLatcdfs=figure;
    subplot(2,1,1);
    [rftF,rftX,rftO,rftUp]=ecdf(lickLatRFT(:),'Bounds','on');
    [fulltF,fulltX,fulltO,fulltUp]=ecdf(lickLatFullT(:),'Bounds','on');
    [ttF,ttX,ttO,ttUp]=ecdf(lickLatToneT(:),'Bounds','on');
    allLickLatCFullT=[lickLatCFullT;lickLatCT];
    [ctF,ctX,ctO,ctUp]=ecdf(allLickLatCFullT(:),'Bounds','on');
    plot(rftX,rftF,'Color',reinfcolor,'LineWidth',2);hold on;
    plot(rftX,rftO,':','Color',reinfcolor,'LineWidth',1);plot(rftX,rftUp,':','Color',reinfcolor,'LineWidth',1);
    plot(fulltX,fulltF,'Color',darkOptoColor,'LineWidth',2);
    plot(fulltX,fulltO,':','Color',darkOptoColor,'LineWidth',1);plot(fulltX,fulltUp,':','Color',darkOptoColor,'LineWidth',1);
    plot(ttX,ttF,'Color',optocolor,'LineWidth',2);
    plot(ttX,ttO,':','Color',optocolor,'LineWidth',1); plot(ttX,ttUp,':','Color',optocolor,'LineWidth',1);
    plot(ctX,ctF,'Color',lighterOptoColor,'LineWidth',2);
    plot(ctX,ctO,':','Color',lighterOptoColor,'LineWidth',1);plot(ctX,ctUp,':','Color',lighterOptoColor,'LineWidth',1);
    
    legend('Light off','','','Full','','','Stimulus','','','Choice','','','Location','Best');xlim([0 4]);
    ylabel('p(first lick)');xlabel('Time of first S+ lick (s)');box off;
    summaryLickLatcdfs.Position(3:4)=[250 500];
    saveas(gcf,['allmice_T_MGB_HitFull_LickLatSum']);
    saveas(gcf,['allmice_T_MGB_HitFull_LickLatSum.png']);
    
    
    subplot(2,1,2);hold off;
    [rffF,rffX,rffO,rffUp]=ecdf(lickLatRFF(:),'Bounds','on');
    [fullfF,fullfX,fullfO,fullfUp]=ecdf(lickLatFullF(:),'Bounds','on');
    [tfF,tfX,tfO,tfUp]=ecdf(lickLatToneF(:),'Bounds','on');
    allLickLatCFullF=[lickLatCFullF;lickLatCF];
    [cfF,cfX,cfO,cfUp]=ecdf(allLickLatCFullF(:),'Bounds','on');
    plot(rffX,rffF,'Color',reinfcolor,'LineWidth',2);hold on;
    plot(rffX,rffO,':','Color',reinfcolor,'LineWidth',1);
    plot(rffX,rffUp,':','Color',reinfcolor,'LineWidth',1);
    plot(fullfX,fullfF,'Color',darkOptoColor,'LineWidth',2);
    plot(fullfX,fullfO,':','Color',darkOptoColor,'LineWidth',1);
    plot(fullfX,fullfUp,':','Color',darkOptoColor,'LineWidth',1);
    plot(tfX,tfF,'Color',optocolor,'LineWidth',2);
    plot(tfX,tfO,':','Color',optocolor,'LineWidth',1); plot(tfX,tfUp,':','Color',optocolor,'LineWidth',1);
    plot(cfX,cfF,'Color',lighterOptoColor,'LineWidth',2);
    plot(cfX,cfO,':','Color',lighterOptoColor,'LineWidth',1);
    plot(cfX,cfUp,':','Color',lighterOptoColor,'LineWidth',1);
    legend('Light off','','','Full','','','Stimulus','','','Choice','','','Location','Best');xlim([0 4]);
    ylabel('p(first lick)');xlabel('Time of first S+ lick (s)');box off;
    summaryLickLatcdfs.Position(3:4)=[250 500];
    saveas(gcf,['allmice_T_MGB_ecdf_LickLatSum']);
    saveas(gcf,['allmice_T_MGB_ecdf_LickLatSum.pdf']);
    
    
    % Lick latency PSTHs
    lickLatHistsAll=figure; % separated by pairs of same-day condition 
    subplot(1,4,1);hold on;title('Hit Trials');
    plot(bins,RFTHist,'Color',reinfcolor,'LineWidth',2);
    plot(bins,FullTHist,'Color',(optocolor-0.3),'LineWidth',2);xlabel('Time of first lick(s)');
    plot(bins,CFTHist,'Color',lightOptoColor,'LineWidth',2);xlabel('Time of first lick(s)');
    legend('light off','full','choice');ylabel('Proportion of licks');
    ylim([0 0.5]);xlim([-0.5 4]);
    subplot(1,4,2);hold on;title('Hit Trials');
    plot(bins,RTTHist,'Color',reinfcolor,'LineWidth',2);
    plot(bins,ToneTHist,'Color',(optocolor-0.3),'LineWidth',2);xlabel('Time of first lick(s)');
    plot(bins,CTHist,'Color',lightOptoColor,'LineWidth',2);xlabel('Time of first lick(s)');
    legend('light off','stimulus','choice');
    ylim([0 0.5]);xlim([-0.5 4]);
    subplot(1,4,3);hold on;title('FA Trials');
    plot(bins,RFFHist,'Color',reinfcolor,'LineWidth',2);
    plot(bins,FullFHist,'Color',(optocolor-0.3),'LineWidth',2);xlabel('Time of first lick(s)');
    plot(bins,CFFHist,'Color',lightOptoColor,'LineWidth',2);xlabel('Time of first lick(s)');
    legend('light off','full','choice');
    ylim([0 0.5]);xlim([-0.5 4]);
    subplot(1,4,4);hold on;title('FA Trials');
    plot(bins,RTFHist,'Color',reinfcolor,'LineWidth',2);
    plot(bins,ToneFHist,'Color',(optocolor-0.3),'LineWidth',2);xlabel('Time of first lick(s)');
    plot(bins,CFHist,'Color',lightOptoColor,'LineWidth',2);xlabel('Time of first lick(s)');
    ylim([0 0.5]);xlim([-0.5 4]);
    lickLatHistsAll.Position(3:4)=[800 250];
    
    
    lickLatHistsAllTogether=figure; % separated by pairs of same-day condition 
    subplot(1,2,1);hold on;title('Hit Trials');
    plot(bins,(RFTHist+RTTHist)/2,'Color',reinfcolor,'LineWidth',2);
    plot(bins,FullTHist,'Color',(optocolor-0.4),'LineWidth',2);xlabel('Time of first lick(s)');
    plot(bins,ToneTHist,'Color',(optocolor-0.1),'LineWidth',2);
    plot(bins,(CFTHist+CTHist)/2,'Color',lighterOptoColor,'LineWidth',2);xlabel('Time of first lick(s)'); ylabel('Proportion of licks');
%     plot(bins,CTHist,'Color',lightOptoColor,'LineWidth',2);xlabel('Time of first lick(s)');
    legend('light off','full','stimulus','choice');
    ylim([0 0.6]);xlim([-0.5 4]);
%     aucTargetN=trapz(bins,(RFTHist+RTTHist)/2);
%     aucTargetFull=trapz(bins,FullTHist);
%     aucTargetTone=trapz(bins,ToneTHist);
%     aucTargetChoice=trapz(bins,(CFTHist+CTHist)/2);
    
    %%%%%%%% need to calcualte AUC by doing it per animal and then
    %%%%%%%% averaging to get error bars
%     subplot(2,2,2); % add bar plot here quantifying AUC
%     hitBar=bar([aucTargetN; aucTargetFull; aucTargetTone; aucTargetChoice],'LineWidth',1.5);hold on;
%     hitBar(1).FaceColor='none';hitBar(1).EdgeColor='flat'; 
%     hitBar.CData(1,:)=reinfcolor;hitBar.CData(2:4,:)=optocolor;
    
    subplot(1,2,2);hold on;title('FA Trials');
    plot(bins,(RFFHist+RTFHist)/2,'Color',reinfcolor,'LineWidth',2);
    plot(bins,FullFHist,'Color',(optocolor-0.4),'LineWidth',2);xlabel('Time of first lick(s)');
    plot(bins,ToneFHist,'Color',(optocolor-0.1),'LineWidth',2);xlabel('Time of first lick(s)');
    plot(bins,(CFHist+CFFHist)/2,'Color',lighterOptoColor,'LineWidth',2);xlabel('Time of first lick(s)');
    ylim([0 0.5]);xlim([-0.5 4]);
    legend('light off','full','stimulus','choice');
    ylim([0 0.6]);xlim([-0.5 4]);
    lickLatHistsAllTogether.Position(3:4)=[300 200];
    
    trialData{1,1}='Condition, Condition Day';trialData{1,2}='Target';trialData{1,4}='Foil';
    trialData{2,1}='No Light, Full Trial Day';trialData{2,2}=allMiceRFTarget;
    trialData{3,1}='Full Trial, Full Trial Day';trialData{3,2}=allMiceFullTarget;
    trialData{2,4}=allMiceRFFoil;trialData{3,4}=allMiceFullFoil;
    trialData{4,1}='Choice, Full Trial Day';trialData{4,2}=allMiceCFTarget;trialData{4,4}=allMiceCFFoil;
    trialData{5,1}='No Light, Tone Day';trialData{5,2}=allMiceRTTarget;trialData{5,4}=allMiceRTFoil;
    trialData{6,1}='Tone, Tone Day';trialData{6,2}=allMiceToneTarget;trialData{6,4}=allMiceToneFoil;
    trialData{7,1}='Choice, Tone Day';trialData{7,2}=allMiceCTTarget;trialData{7,4}=allMiceCTFoil;
        trialData{1,3}='Target Idx';
    trialData{2,3}=idxAllMice{7,3};
    trialData{3,3}=idxAllMice{9,3};
    trialData{4,3}=idxAllMice{11,3};
    trialData{5,3}=idxAllMice{1,3};
    trialData{6,3}=idxAllMice{3,3};
    trialData{7,3}=idxAllMice{5,3};
    trialData{1,5}='Foil Idx';
    trialData{2,5}=idxAllMice{8,3};
    trialData{3,5}=idxAllMice{10,3};
    trialData{4,5}=idxAllMice{12,3};
    trialData{5,5}=idxAllMice{2,3};
    trialData{6,5}=idxAllMice{4,3};
    trialData{7,5}=idxAllMice{6,3};
    
    save('byTrialData','trialData');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% all animals, rolling average %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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

        fullTargetIdx=movmean(mean(allMiceDelay1TargetIdx,1),windowsize)';
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