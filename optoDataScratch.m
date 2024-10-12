% load the data
cd('O:\sjk\Behavior\OptoMGBIC_1')
load('summaryData.mat');
cohort1=optomeanMat;
cd('O:\sjk\Behavior\OptoMGBIC_2')
load('summaryData.mat');
cohort2=optomeanMat;
cd('O:\sjk\Behavior\OptoMGBIC_3')
load('summaryData.mat');
cohort3=optomeanMat;
cd('O:\sjk\Behavior\MGBIC_4')
load('summaryData.mat');
cohort4=optomeanMat;

% cd('O:\sjk\Behavior\MGBIC_5')
% load('summaryData.mat');
% cohort5=optomeanMat;


% allCohorts=cohort1;
% dim=size(allCohorts);
% dimToAdd=size(cohort2);
% allCohorts(dim+1:dim+dimToAdd(1)-1,1:dim(2))=cohort2(2:dimToAdd(1),1:dim(2));
% dimToAdd=size(cohort3);
% dim=size(allCohorts);
% allCohorts(dim(1)+1:dim(1)+dimToAdd(1)-1,1:dim(2))=cohort3(2:dimToAdd(1),1:dim(2));
% dimToAdd=size(cohort4);
% dim=size(allCohorts);
% allCohorts(dim(1)+1:dim(1)+dimToAdd(1)-1,1:dim(2))=cohort4(2:dimToAdd(1),1:dim(2));
% save('allOptoCohortData.mat','allCohorts','cohort1','cohort2','cohort3','cohort4','cohort5');

%% here is to load the data from the most recent summary file for cohorts 1-4
% load the summary data
cd('C:\Users\sjkim1\Desktop\OptoData')
load('allOptoCohortData.mat');

% add cohort 5
cd('C:\Users\sjkim1\Desktop\OptoData\MGBIC_5')
load('summaryData.mat');
cohort5=optomeanMat;
cd('C:\Users\sjkim1\Desktop\OptoData');
dimToAdd=size(cohort5);
dim=size(allCohorts);
allCohorts(dim(1)+1:dim(1)+dimToAdd(1)-1,1:dim(2))=cohort5(2:dimToAdd(1),1:dim(2));
save('allOptoCohortData.mat','allCohorts','cohort5');
%%
% compute cohen's d
% medium effect of 0.5 is visible to the naked eye
% condition=[0 0 1 0 1 0 1 0 1 ...
%     0 2 0 2 0 1 0 2 ...
%     0 2 0 1 0 2 0 1]; % 1 is ctl, 2 is test
condition=[0 0 1 0 1 0 1 0 1 ...
    0 2 0 1 0 1 0 2 ... % 164 is a test btu coding as a CTL since little effect/high baseline FA rate
    0 2 0 1 0 1]; % 1 is ctl, 2 is test
testIdx=find(condition==2);
SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
ii=1;
    
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

binRateSize=5;ee=1;
cohortNum = [1 2 3]; % for test only 
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

%% cohens D test    
    
cohensDTest={}; % this will give you effect sizes we see in the pilot data
% using the mean of each session
cohensDTest{1,1}={'Animal Name'};
cohensDTest{1,2}={'MGB Full Opto Hit'};
cohensDTest{1,3}={'MGB Full Opto FA'};
cohensDTest{1,4}={'MGB Tone Opto Hit'};
cohensDTest{1,5}={'MGB Tone Opto FA'};
cohensDTest{1,6}={'MGB Choice Opto Hit'};
cohensDTest{1,7}={'MGB Choice Opto FA'};
pp=2;

% correct for high and low values
% newBinRates is uncorrected, binRates is corrected

% newBinRates=binRates;
uu=2;yy=3;
for uu=2:size(binRates,1)
    for yy=3:size(binRates,2)
        dimSub=length(binRates{uu,yy});
        if find(cell2mat(binRates(uu,yy))==1)
            idx=find(cell2mat(binRates(uu,yy))==1);
            binRates{uu,yy}(idx)=(dimSub-1)/dimSub;
        elseif find(cell2mat(binRates(uu,yy))==0)
            idx=find(cell2mat(binRates(uu,yy))==0);
            binRates{uu,yy}(idx)=1/dimSub;
        else
        end
    end
end

for tt=2:size(binRates,1)
    cohensDTest{pp,1}=binRates{tt,1};
    cohensDTest{pp,2}=computeCohen_d(binRates{tt,3}(1:(length(binRates{tt,5}))),...
        binRates{tt,5},'paired'); % Full Opto H
    mgbRHitStd=std(binRates{tt,3}(1:(length(binRates{tt,5}))));
    cohensDTest{pp,3}=computeCohen_d(binRates{tt,4}(1:(length(binRates{tt,6}))),...
        binRates{tt,6},'paired'); % Full Opto FA
    mgbRFaStd=std(binRates{tt,4}(1:(length(binRates{tt,6}))));
    cohensDTest{pp,4}=computeCohen_d(binRates{tt,3}(1:(length(binRates{tt,7}))),...
        binRates{tt,7},'paired'); % Tone Opto H
    cohensDTest{pp,5}=computeCohen_d(binRates{tt,4}(1:(length(binRates{tt,8}))),...
        binRates{tt,8},'paired'); % Tone Opto FA
    cohensDTest{pp,6}=computeCohen_d(binRates{tt,3}(1:(length(binRates{tt,9}))),...
        binRates{tt,9},'paired'); % Choice Opto H
    cohensDTest{pp,7}=computeCohen_d(binRates{tt,4}(1:(length(binRates{tt,10}))),...
        binRates{tt,10},'paired'); % Choice Opto FA
    pp=pp+1;
end
%%    % how many samples do i need for each value??? % Inputs for sample size est.
nn=1:20; count=-1;
tt=2;
figure(tt);
for tt=2:size(binRates,1)
     %hits first
    subplot(3,3,tt+count);
    nout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
        mean(binRates{tt,5}),0.80); % % (type of test, [mean1 stdev], mean2, power)
    pwrout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
        mean(binRates{tt,5}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),[],nn);

    nout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),[],nn);
    plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');ylim([0.65 1]);
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    if tt==4
        xlim([0 10]);
    else
    end
    title([char(binRates{tt,1}) 'MGB Hit Full Power vs. Sample size']);
    xlabel('sample size (bins of 5 trials)');ylabel('power');
    
    nout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),[],nn);
    count=count+1; subplot(3,3,tt+count);plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');xlim([0 10]);
    title('Hit Tone Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');

    nout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),[],nn); 
    count=count+1;subplot(3,3,tt+count); plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');xlim([0 10]);
    title('Hit Choice Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');
end

nn=1:55; count=-1;
figure(tt+1);
for tt=2:size(binRates,1)
    subplot(3,3,tt+count);
    nout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),0.80);
    pwrout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})))) mgbRFaStd], ...
    mean(binRates{tt,6}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),[],nn); 
    plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');ylim([0.65 1]);
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    title([char(binRates{tt,1}) 'MGB FA Full Power vs. Sample size']);
    xlabel('sample size (bins of 5 trials)');ylabel('power');
    
    nout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),[],nn);
    count=count+1; subplot(3,3,tt+count);plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');
    title('MGB FA Tone Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');
    
    nout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),[],nn); 
    count=count+1;subplot(3,3,tt+count); 
    plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');
    title('MGB FA Choice Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');
end


%%

% for a given animal, take effect sizes from cohensDTest variable and then
% calculate the sample sizes needed. 
for uu=2:size(cohensDTest,1)-2
    effect_sizes = [cohensDTest{uu,2}, cohensDTest{uu,3},cohensDTest{uu,4},cohensDTest{uu,5},...
        cohensDTest{uu,6},cohensDTest{uu,7}];  
    effect_str   = {'MGB Full Hit', 'MGB Full FA','MGB Tone Hit','MGB Tone FA',...
        'MGB Choice Hit','MGB Choice FA'};
    alpha        = 0.05; % 0.05, 0.01, 0.005, 0.001                   
    stat_power   = 0.85; % 0.80, 0.85, 0.90                    

    % Sample size estimates
    sample_sizes = arrayfun(@(effect) calculateSampleSize(effect, alpha, stat_power), effect_sizes);
%     sample_sizes = arrayfun(@(effect) sampsizepwr('t2', [effect 1], alpha, stat_power), effect_sizes);

    for i = 1:length(effect_sizes)
        fprintf('For shift between %s conditions an effect size of %.4f requires a sample size of %d\n', effect_str{i}, effect_sizes(i), ceil(sample_sizes(i)));
    end

end

%%
% catch trials
% coded as conffa (catch light on false alarm
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
%%

condition=[0 0 1 0 1 0 1 0 1 ...
    0 2 0 1 0 1 0 2 ... % 164 is a test btu coding as a CTL since little effect/high baseline FA rate
    0 2 0 1 0 1 0 1 0 1 0 1]; % 1 is ctl, 2 is test
testIdx=find(condition==2);

allDataTestsOnly=allCohorts(1,:);
allDataTestsOnly(2:4,:)=allCohorts(testIdx,:);

ctlIdx=[21,25,27,31];
allDataCtlOnly=allCohorts(1,:);
allDataCtlOnly(2:5,:)=allCohorts(ctlIdx,:);

reinfcolor= [0.4,0.4,0.4];
optocolor=[102/255 178/255 255/255];
    
%%
% make plto to compare percentage correct when light is on vs. off
% Compute percent correct for each session
allDataTestsOnly{1,27}='RPC MGB Full Trial';
allDataTestsOnly{1,28}='OPC MGB Full Trial';
allDataTestsOnly{1,29}='RPC MGB Tone Trial';
allDataTestsOnly{1,30}='OPC MGB Tone Trial';
allDataTestsOnly{1,31}='RPC MGB Choice Trial';
allDataTestsOnly{1,32}='OPC MGB Choice Trial';

for jj=2:4
    eeFig=figure(jj);hold on;
    rhit=allDataTestsOnly{jj,10}; % full trial MGB
    rfa=allDataTestsOnly{jj,11};
    ohit = allDataTestsOnly{jj,12};
    ofa = allDataTestsOnly{jj,13};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(1,3,1)
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
    subplot(1,3,2)
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
    
    subplot(1,3,3)
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
    eeFig.Position(3:4)=[725 275];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_PercentCorrect_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_PercentCorrect_Opto.png']);
end


allDataCtlOnly{1,27}='RPC MGB Full Trial';
allDataCtlOnly{1,28}='OPC MGB Full Trial';
allDataCtlOnly{1,29}='RPC MGB Tone Trial';
allDataCtlOnly{1,30}='OPC MGB Tone Trial';
allDataCtlOnly{1,31}='RPC MGB Choice Trial';
allDataCtlOnly{1,32}='OPC MGB Choice Trial';

for jj=2:5
    eeFig=figure(jj+10);hold on;
    rhit=allDataCtlOnly{jj,10}; % full trial MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,12};
    ofa = allDataCtlOnly{jj,13};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    rpc = (rhit+(1-rfa))/2*100; 
    opc = (ohit+(1-ofa))/2*100; 
    subplot(1,3,1)
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
    subplot(1,3,2)
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
    
    subplot(1,3,3)
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
    eeFig.Position(3:4)=[725 275];
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_PercentCorrect_Opto']);
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_PercentCorrect_Opto.png']);
end

%%
% make plto to compare hit and fa rate with light on vs light off
% for jj=2:4
%     eeFig=figure(jj+3);hold on;
%     rhit=allDataTestsOnly{jj,10}; % full trial MGB
%     rfa=allDataTestsOnly{jj,11};
%     ohit = allDataTestsOnly{jj,12};
%     ofa = allDataTestsOnly{jj,13};
%     rhit(isnan(ohit))=nan;
%     rfa(isnan(ofa))=nan;
%     subplot(1,3,1)
%     eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
%     eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
%         rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
%     scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
%         ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
%     scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
%         rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
%     scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
%         ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
%     [h,pHit,ci,stats] = ttest2(rhit,ohit);
%     [h,pFA,ci,stats] = ttest2(rfa,ofa);
%     sigstar({[1,2],[3,4]}, [pHit pFA])
%     ylabel('rate');
%     title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
%     xticklabels({'hit', 'hit','fa','fa'});
% 
%     rhit=allDataTestsOnly{jj,10}; % tone MGB
%     rfa=allDataTestsOnly{jj,11};
%     ohit = allDataTestsOnly{jj,14};
%     ofa = allDataTestsOnly{jj,15};
%     rhit(isnan(ohit))=nan;
%     rfa(isnan(ofa))=nan;
%     subplot(1,3,2)
%     eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
%     eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
%         rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
%     scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
%         ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
%     scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
%         rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
%     scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
%         ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
%     [h,pHit,ci,stats] = ttest2(rhit,ohit);
%     [h,pFA,ci,stats] = ttest2(rfa,ofa);
%     sigstar({[1,2],[3,4]}, [pHit pFA])
%     ylabel('rate');
%     title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
%     xticklabels({'hit', 'hit','fa','fa'});
% 
%     rhit=allDataTestsOnly{jj,10}; % choice MGB
%     rfa=allDataTestsOnly{jj,11};
%     ohit = allDataTestsOnly{jj,16};
%     ofa = allDataTestsOnly{jj,17};
%     rhit(isnan(ohit))=nan;
%     rfa(isnan(ofa))=nan;
%     rpc = (rhit+(1-rfa))/2*100; 
%     opc = (ohit+(1-ofa))/2*100; 
%     subplot(1,3,3)
%     eee=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
%     eee(1).FaceColor='flat'; eee(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     scatter(repmat(eee(1).XEndPoints(1),size(rhit,1),1), ...
%         rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
%     scatter(repmat(eee(1).XEndPoints(2),size(ohit,1),2), ...
%         ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
%     scatter(repmat(eee(1).XEndPoints(3),size(rfa,1),1), ...
%         rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
%     scatter(repmat(eee(1).XEndPoints(4),size(ofa,1),2), ...
%         ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
%     [h,pHit,ci,stats] = ttest2(rhit,ohit);
%     [h,pFA,ci,stats] = ttest2(rfa,ofa);
%     sigstar({[1,2],[3,4]}, [pHit pFA])
%     ylabel('rate');
%     title([allDataTestsOnly{jj,1} ' MGB Choice Inactivation']);
%     xticklabels({'hit', 'hit','fa','fa'});
%     eeFig.Position(3:4)=[725 275];
%     saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_HitFARate_Opto']);
%     saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_HitFARate_Opto.png']);
% end

%CONTROLS
for jj=2:5
    eeFig=figure(jj+3);hold on;
    rhit=allDataCtlOnly{jj,10}; % full trial MGB
    rfa=allDataCtlOnly{jj,11};
    ohit = allDataCtlOnly{jj,12};
    ofa = allDataCtlOnly{jj,13};
    rhit(isnan(ohit))=nan;
    rfa(isnan(ofa))=nan;
    subplot(1,3,1)
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
    subplot(1,3,2)
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
    subplot(1,3,3)
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
    eeFig.Position(3:4)=[725 275];
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_HitFARate_Opto']);
    saveas(gcf,[char(allDataCtlOnly{jj,1}) '_C_MGB_HitFARate_Opto.png']);
end


%% all animal summary analysis
% group by animal
% test
clear rpc opc
for jj=2:size(allDataTestsOnly,1)
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,27});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,28});
end
wwFig=figure(10);
subplot(1,3,1)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;

% 
clear rpc opc
for jj=2:size(allDataTestsOnly,1)
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,29});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,30});
end
subplot(1,3,2)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
for jj=2:size(allDataTestsOnly,1)
    rpc(jj-1)=nanmean(allDataTestsOnly{jj,31});
    opc(jj-1)=nanmean(allDataTestsOnly{jj,32});
end
subplot(1,3,3)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});
wwFig.Position(3:4)=[725 275];
% saveas(gcf,['ByAnimal_T_MGB_PercentCorrect_Opto']);
% saveas(gcf,['ByAnimal_T_MGB_PercentCorrect_Opto.png']);

% Controls
clear rpc opc
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,27});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,28});
end
wwFig=figure(10);
subplot(1,3,1)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataCtlOnly{jj,27}=rpc;
allDataCtlOnly{jj,28}=opc;


clear rpc opc
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,29});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,30});
end
subplot(1,3,2)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
for jj=2:size(allDataCtlOnly,1)
    rpc(jj-1)=nanmean(allDataCtlOnly{jj,31});
    opc(jj-1)=nanmean(allDataCtlOnly{jj,32});
end
subplot(1,3,3)
qqq=bar([mean(rpc) mean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});
wwFig.Position(3:4)=[725 275];
saveas(gcf,['ByAnimal_C_MGB_PercentCorrect_Opto']);
saveas(gcf,['ByAnimal_C_MGB_PercentCorrect_Opto.png']);

%% now group by session, percent correct for Test
reinfcolor1= [0.2,0.2,0.2];
reinfcolor2= [0.45,0.45,0.45];
reinfcolor3= [0.7,0.7,0.7];
optocolor1=[0/255 89/255 178/255];
optocolor2=[65/255 161/255 255/255];
optocolor3=[181/255 216/255 255/255];

clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,27});
    opc=cat(1,opc,allDataTestsOnly{jj,28});
end
wwFig=figure(10);
subplot(1,3,1)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataTestsOnly{jj,27}=rpc;
allDataTestsOnly{jj,28}=opc;


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,29});
    opc=cat(1,opc,allDataTestsOnly{jj,30});
end
subplot(1,3,2)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:13),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,31});
    opc=cat(1,opc,allDataTestsOnly{jj,32});
end
subplot(1,3,3)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor2);
scatter(repmat(qqq(1).XEndPoints(1),4,1), ...
    rpc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor3);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor1);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(6:9),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor2);
scatter(repmat(qqq(1).XEndPoints(2),4,1), ...
    opc(10:16),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor3);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});
wwFig.Position(3:4)=[725 275];
saveas(gcf,['BySess_T_MGB_PercentCorrect_Opto']);
saveas(gcf,['BySess_T_MGB_PercentCorrect_Opto.png']);

%% controls by session percent correect
allDataCtlOnly{5,27}=allDataCtlOnly{5,27}';
allDataCtlOnly{5,28}=allDataCtlOnly{5,28}';

clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,27});
    opc=cat(1,opc,allDataCtlOnly{jj,28});
end
wwFig=figure(10);
subplot(1,3,1)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
ylabel('percent correct');
title(['MGB Full Trial Inactivation']);
xticklabels({'light off', 'light on'});
allDataCtlOnly{jj,27}=rpc;
allDataCtlOnly{jj,28}=opc;


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,29});
    opc=cat(1,opc,allDataCtlOnly{jj,30});
end
subplot(1,3,2)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Tone Inactivation']);
xticklabels({'light off', 'light on'});


clear rpc opc
rpc=NaN;opc=NaN;
for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,31});
    opc=cat(1,opc,allDataCtlOnly{jj,32});
end
subplot(1,3,3)
qqq=bar([nanmean(rpc) nanmean(opc)]); hold on;
qqq(1).FaceColor='flat'; qqq(1).CData=[reinfcolor;optocolor];hold on;
scatter(repmat(qqq(1).XEndPoints(1),size(rpc,1),1), ...
    rpc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qqq(1).XEndPoints(2),size(opc,1),1), ...
    opc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(rpc,opc);
sigstar({[1,2]}, p)
title(['MGB Choice Inactivation']);
xticklabels({'light off', 'light on'});
wwFig.Position(3:4)=[725 275];
saveas(gcf,['BySess_C_MGB_PercentCorrect_Opto']);
saveas(gcf,['BySess_C_MGB_PercentCorrect_Opto.png']);

 %% now do by animal for hit and FA for Test animals

ppFig=figure(15);

for jj=2:size(allDataTestsOnly,1)
    rpc=cat(1,rpc,allDataTestsOnly{jj,27});
    opc=cat(1,opc,allDataTestsOnly{jj,28});
end

clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,10}; % full MGB
    rfatemp=allDataTestsOnly{jj,11};
    ohittemp=allDataTestsOnly{jj,12};
    ofatemp=allDataTestsOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(1,3,1)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,2),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,10}; % tone MGB
    rfatemp=allDataTestsOnly{jj,11};
    ohittemp = allDataTestsOnly{jj,14};
    ofatemp = allDataTestsOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(1,3,2)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,2),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,10}; % choice MGB
    rfatemp=allDataTestsOnly{jj,11};
    ohittemp = allDataTestsOnly{jj,16};
    ofatemp = allDataTestsOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(1,3,3)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,2),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});
ppFig.Position(3:4)=[725 275];
saveas(gcf,['By Animal T_MGB_HitFARate_Opto']);
saveas(gcf,['By Animal T_MGB_HitFARate_Opto.png']);


%% CONTROL by animal hit and false alarm rate

allDataCtlOnly{5,27}=allDataCtlOnly{5,27}';
allDataCtlOnly{5,28}=allDataCtlOnly{5,28}';
ppFig=figure(15);

for jj=2:size(allDataCtlOnly,1)
    rpc=cat(1,rpc,allDataCtlOnly{jj,27});
    opc=cat(1,opc,allDataCtlOnly{jj,28});
end

clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; % full MGB
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,12};
    ofatemp=allDataCtlOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(1,3,1)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,2),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; % tone MGB
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp = allDataCtlOnly{jj,14};
    ofatemp = allDataCtlOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(1,3,2)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,2),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; % choice MGB
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp = allDataCtlOnly{jj,16};
    ofatemp = allDataCtlOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit(jj-1)=nanmean(rhittemp);
    rfa(jj-1)=nanmean(rfatemp);
    ohit(jj-1)=nanmean(ohittemp);
    ofa(jj-1) = nanmean(ofatemp);
end
subplot(1,3,3)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,2),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,2),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,2),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,2),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});
ppFig.Position(3:4)=[725 275];
saveas(gcf,['By Animal C_MGB_HitFARate_Opto']);
saveas(gcf,['By Animal C_MGB_HitFARate_Opto.png']);

%% By session H and FA
reinfcolor1= [0.2,0.2,0.2];
reinfcolor2= [0.45,0.45,0.45];
reinfcolor3= [0.7,0.7,0.7];
optocolor1=[0/255 89/255 178/255];
optocolor2=[65/255 161/255 255/255];
optocolor3=[181/255 216/255 255/255];

ppFig=figure(12);
clear rhit ohit rfa ofa
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,10}; 
    rfatemp=allDataTestsOnly{jj,11};
    ohittemp=allDataTestsOnly{jj,12};
    ofatemp=allDataTestsOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(1,3,1)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor1);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit(1:4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,10}; 
    rfatemp=allDataTestsOnly{jj,11};
    ohittemp=allDataTestsOnly{jj,14};
    ofatemp=allDataTestsOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(1,3,2)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataTestsOnly,1)
    rhittemp=allDataTestsOnly{jj,10}; 
    rfatemp=allDataTestsOnly{jj,11};
    ohittemp=allDataTestsOnly{jj,16};
    ofatemp=allDataTestsOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(1,3,3)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});
ppFig.Position(3:4)=[725 275];
saveas(gcf,['By Sess T_MGB_HitFARate_Opto']);
saveas(gcf,['By Sess T_MGB_HitFARate_Opto.png']);

%% by session, hit and false alarm rate Controls

ppFig=figure(12);
clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; 
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,12};
    ofatemp=allDataCtlOnly{jj,13};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(1,3,1)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['By Animal MGB Full Trial Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; 
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,14};
    ofatemp=allDataCtlOnly{jj,15};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(1,3,2)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title(['MGB Tone Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});

clear rhit ohit rfa ofa ohittemp rhittemp fahittemp ofatemp
rhit=NaN;ohit=NaN;rfa=NaN;ofa=NaN;
for jj=2:size(allDataCtlOnly,1)
    rhittemp=allDataCtlOnly{jj,10}; 
    rfatemp=allDataCtlOnly{jj,11};
    ohittemp=allDataCtlOnly{jj,16};
    ofatemp=allDataCtlOnly{jj,17};
    rhittemp(isnan(ohit))=nan;
    rhittemp(rhittemp==0)=nan;
    rfatemp(isnan(ofa))=nan;
    rhit=cat(1,rhit,rhittemp);
    rfa=cat(1,rfa,rfatemp);
    ohit=cat(1,ohit,ohittemp);
    ofa = cat(1,ofa,ofatemp);
end
subplot(1,3,3)
ppp=bar([nanmean(rhit) nanmean(ohit) nanmean(rfa) nanmean(ofa)]); hold on;
ppp(1).FaceColor='flat'; ppp(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
scatter(repmat(ppp(1).XEndPoints(1),size(rhit,1),1), ...
    rhit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(2),size(ohit,1),2), ...
    ohit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(ppp(1).XEndPoints(3),size(rfa,1),1), ...
    rfa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(ppp(1).XEndPoints(4),size(ofa,1),2), ...
    ofa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,pHit,ci,stats] = ttest2(rhit,ohit);
[h,pFA,ci,stats] = ttest2(rfa,ofa);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('rate');
title([' MGB Choice Inactivation']);
xticklabels({'hit', 'hit','fa','fa'});
ppFig.Position(3:4)=[725 275];
saveas(gcf,['By Sess C_MGB_HitFARate_Opto']);
saveas(gcf,['By Sess C_MGB_HitFARate_Opto.png']);

%% now per animal statistical analysis as a 2way anova
% no light vs. light all inactivation condition
testFullTrial=NaN(18,2); 
% Full Trial, light OFF column 1, light ON column 2
testFullTrial(1:2,1)=allDataTestsOnly{2,27}(~isnan(allDataTestsOnly{2,27}));
testFullTrial(1:2,2)=allDataTestsOnly{2,28}(~isnan(allDataTestsOnly{2,28}));
testFullTrial(3:4,1)=allDataTestsOnly{3,27}(~isnan(allDataTestsOnly{3,27}));
testFullTrial(3:4,2)=allDataTestsOnly{3,28}(~isnan(allDataTestsOnly{3,28}));
fullTrialtemp=allDataTestsOnly{4,27}(~isnan(allDataTestsOnly{4,27}));
testFullTrial(5:6,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,28}(~isnan(allDataTestsOnly{4,28}));
testFullTrial(5:6,2)=fullTrialtemp(1:2);
% tone Trial, light OFF column 1, light ON column 2
testFullTrial(7:8,1)=allDataTestsOnly{2,29}(~isnan(allDataTestsOnly{2,29}));
testFullTrial(7:8,2)=allDataTestsOnly{2,30}(~isnan(allDataTestsOnly{2,30}));
testFullTrial(9:10,1)=allDataTestsOnly{3,29}(~isnan(allDataTestsOnly{3,29}));
testFullTrial(9:10,2)=allDataTestsOnly{3,30}(~isnan(allDataTestsOnly{3,30}));
fullTrialtemp=allDataTestsOnly{4,29}(~isnan(allDataTestsOnly{4,29}));
testFullTrial(11:12,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,30}(~isnan(allDataTestsOnly{4,30}));
testFullTrial(11:12,2)=fullTrialtemp(1:2);
% Choice Trial, light OFF column 1, light ON column 2
fullTrialtemp=allDataTestsOnly{2,31}(~isnan(allDataTestsOnly{2,31}));
testFullTrial(13:14,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{2,32}(~isnan(allDataTestsOnly{2,32}));
testFullTrial(13:14,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{3,31}(~isnan(allDataTestsOnly{3,31}));
testFullTrial(15:16,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{3,32}(~isnan(allDataTestsOnly{3,32}));
testFullTrial(15:16,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,31}(~isnan(allDataTestsOnly{4,31}));
testFullTrial(17:18,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,32}(~isnan(allDataTestsOnly{4,32}));
testFullTrial(17:18,2)=fullTrialtemp(1:2);

[aov,~,stats]=anova2(testFullTrial,6,FactorNames=["Light on vs. Off","Opto Condition"]);
c1 = multcompare(stats);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2 = multcompare(stats,"Estimate","row");
tbl2 = array2table(c2,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
figure;
light = [repmat("On",18,1);repmat("Off",18,1);];
optoCond = [repmat("Full Trial",6,1);repmat("Tone",6,1);repmat("Choice",6,1)];
factors = {light, optoCond};
aov = anova(factors,testFullTrial(:),FactorNames=["Light","optoCond"],...
    ModelSpecification="interactions");
boxchart(aov,["Light on vs. Off","Opto Condition"])

%% compare differences in percent correct across each light off vs light on condition

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
%%

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

