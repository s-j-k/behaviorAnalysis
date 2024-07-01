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

allCohorts={}; 
allCohorts=cohort1;
dim=size(allCohorts);
dimToAdd=size(cohort2);
allCohorts(dim+1:dim+dimToAdd(1)-1,1:dim(2))=cohort2(2:dimToAdd(1),1:dim(2));
dimToAdd=size(cohort3);
dim=size(allCohorts);
allCohorts(dim(1)+1:dim(1)+dimToAdd(1)-1,1:dim(2))=cohort3(2:dimToAdd(1),1:dim(2));
% need to correct for 1's and 0's in the action rate data; these will make
% the cohens d values later calculate as infinite
save('allOptoCohortData.mat','allCohorts','cohort1','cohort2','cohort3');
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
