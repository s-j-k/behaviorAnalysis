function binRates=makeBinRates(binRatesSize,cohortNum)
% bin the data according to numer of trials specified by binRateSize (5 is
% default) and the number of cohorts (right now, only cohorts 1:3

binRateSize=5;ee=1;
cohortNum = [1 2 3]; % for test only 

% set up the names of the columns 
binRates={};binRates{1,ii}='Animal'; 
binRates{1,ii+1}='Session';
% binRates{1,ii+2}='RHit'; binRates{1,ii+3}='RFA';
% binRates{1,ii+4}='Full Opto Hit'; binRates{1,ii+5}='Full Opto FA';
% binRates{1,ii+6}='Tone Opto Hit'; binRates{1,ii+7}='Tone Opto FA';
% binRates{1,ii+8}='Choice Opto Hit'; binRates{1,ii+9}='Choice Opto FA';
binRates{1,ii+2}= 'MGB R Hit';binRates{1,ii+3}='MGB R FA';
binRates{1,ii+4}='MGB Full O Hit'; binRates{1,ii+5}='MGB Full O FA';
binRates{1,ii+6}='MGB Tone Opto Hit';  binRates{1,ii+7}='MGB Tone Opto FA';
binRates{1,ii+8}='MGB Choice Opto Hit'; binRates{1,ii+9}='MGB Choice Opto FA';
% binRates{1,ii+18}= 'IC R Hit'; binRates{1,ii+19}='IC R FA';
% binRates{1,ii+20}='IC Full Opto Hit';binRates{1,ii+21}='IC Full Opto FA';
% binRates{1,ii+22}='IC Tone Opto Hit';  binRates{1,ii+23}='IC Tone Opto FA';
% binRates{1,ii+24}='IC Choice Opto Hit'; binRates{1,ii+25}='IC Choice Opto FA';
% now use the MAT Variable, which has the full matrix across all days for each animal, 
% to create the binned rates


for ww=1:length(testIdx)
    tt=1;   
    newTempMat=allCohorts{testIdx(ww),26};    % BUT ONLY USE OPTO DAYS
    % need to know which days are the experimental days, and for which
    % condition
    if cohortNum(ww)==1        % for sk163
            mgbDays = [1 1 0 0 1 1 0 0 1 1 1 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 1 1 0 0 1 1 0 0 0 1 1];icDays=logical(icDays);
            expRange=1:8; % Cohort 1
    elseif cohortNum(ww)==2  % sk175
            expRange=22:32; % cohort 2
            mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 0 1 1 0 0 ...
                0 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                1 1 1 0 0 1 1 0 0 1 1 ...
                0 0 0];icDays=logical(icDays);
    elseif cohortNum(ww)==3 %sk176
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ...
                1 0 1 1 0 1 1 0 1 0 0 1];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ... 
                0 1 0 0 1 0 0 1 0 1 1 0];
            icDays=logical(icDays);
            expRange=20:length(mgbDays); 
    end
    mgbDays=mgbDays(expRange);
    icDays=icDays(expRange);
    anyMGB=ismember(newTempMat(:,1),(expRange(mgbDays))); % find MGB opto days....using the first column (day)
    mgbTempMat=newTempMat(anyMGB,:);

    target=mgbTempMat((mgbTempMat(:,CTXT)==2 & mgbTempMat(:,TONE)==1),:);
    otarget=mgbTempMat((mgbTempMat(:,CTXT)==1 & mgbTempMat(:,TONE)==1),:);    
    ottarget=mgbTempMat((mgbTempMat(:,CTXT)==5 & mgbTempMat(:,TONE)==1),:);
    octarget=mgbTempMat((mgbTempMat(:,CTXT)==6 & mgbTempMat(:,TONE)==1),:);
    foil=mgbTempMat((mgbTempMat(:,CTXT)==2 & mgbTempMat(:,TONE)==2),:);
    ofoil=mgbTempMat((mgbTempMat(:,CTXT)==1 & mgbTempMat(:,TONE)==2),:);
    ocfoil=mgbTempMat((mgbTempMat(:,CTXT)==6 & mgbTempMat(:,TONE)==2),:);
    otfoil=mgbTempMat((mgbTempMat(:,CTXT)==5 & mgbTempMat(:,TONE)==2),:);
    numBins=size(target,1)/binRateSize; % bins need to fit the dimensions of how many total trials of that type tehre are
    bincounter=1; obin=1;otbin=1;ocbin=1;
    binSess={};
    binRange=binRateSize;obinRange=binRateSize;otbinRange=binRateSize;ocbinRange=binRateSize;
    yy=1;
    for yy=1:numBins 
        binRhit(1,yy)=sum((target((bincounter:binRange),OUTCOME)==1 & target((bincounter:binRange),CTXT)==2)) / ...
            sum(target((bincounter:binRange),CTXT)==2); % hit in reinforced light off
        binRfa(1,yy) = sum((foil(bincounter:binRange,OUTCOME)==3 & foil(bincounter:binRange,TONE)==2)) / ...
            sum(foil(bincounter:binRange,TONE)==2);
        binSess{1,yy}=[target(bincounter,SESS) target(binRange,SESS)];
        bincounter=bincounter+binRateSize;
        binRange=binRange+binRateSize;
    end
    numBins=size(otarget,1)/binRateSize;yy=1;
    for yy=1:numBins 
        binOhit(1,yy) =  sum((otarget((obin:obinRange),OUTCOME)==1 & otarget((obin:obinRange),CTXT)==1)) / ...
            sum(otarget((obin:obinRange),CTXT)==1); % hit in opto (i.e. reinforced light on)
        binOfa(1,yy) = sum(ofoil(obin:obinRange,OUTCOME)==3 & ofoil(obin:obinRange,CTXT)==1) / ...
            sum(ofoil(obin:obinRange,CTXT)==1);
        obin=obin+binRateSize;
        obinRange=obinRange+binRateSize;
    end
    numBins=size(ottarget,1)/binRateSize;yy=1;
    for yy=numBins
        binOthit(1,yy) = sum((ottarget(otbin:otbinRange,OUTCOME)==1 & ottarget(otbin:otbinRange,CTXT)==5)) / ...
            sum(ottarget(otbin:otbinRange,CTXT)==5); % hit in opto tone (i.e. reinforced light on)
        binOtfa(1,yy) = sum(otfoil(otbin:otbinRange,OUTCOME)==3 & otfoil(otbin:otbinRange,CTXT)==5) / ...
            sum(otfoil(otbin:otbinRange,CTXT)==5); % hit in opto tone (i.e. reinforced light on)
        otbin=otbin+binRateSize;
        otbinRange=otbinRange+binRateSize;
    end
    numBins=size(octarget,1)/binRateSize;yy=1;
    for yy=1:numBins 
        binOchit(1,yy) = sum(octarget(ocbin:ocbinRange,OUTCOME)==1 & octarget(ocbin:ocbinRange,CTXT)==6) / ...
            sum(octarget(ocbin:ocbinRange,CTXT)==6); % hit in opto choice (i.e. reinforced light on)
        binOcfa(1,yy) = sum(ocfoil(ocbin:ocbinRange,OUTCOME)==3 & ocfoil(ocbin:ocbinRange,CTXT)==6) / ...
            sum(ocfoil(ocbin:ocbinRange,CTXT)==6); % hit in opto tone (i.e. reinforced light on)
        ocbin=ocbin+binRateSize;
        ocbinRange=ocbinRange+binRateSize;
    end
    
    binRates{ww+ee,1}=allCohorts{testIdx(ww),1};
    binRates{ww+ee,2}=binSess;
    binRates{ww+ee,3}=binRhit;
    binRates{ww+ee,4}=binRfa; 
    binRates{ww+ee,5}=binOhit;
    binRates{ww+ee,6}=binOfa;
    binRates{ww+ee,7}=binOthit; % is this reinforced light on + light off? 
    binRates{ww+ee,8}=binOtfa;
    binRates{ww+ee,9}=binOchit;
    binRates{ww+ee,10}=binOcfa;
    
end
