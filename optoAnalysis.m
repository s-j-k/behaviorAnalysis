function optoAnalysis(number_mice,pathsave)


nSubj = length(number_mice);
subjlist = cell(1,nSubj);
for i=1:nSubj
    if number_mice(i)-100<0
        subjlist{1,i} = ['SK' num2str(number_mice(i))];
    else
        subjlist{1,i} = ['SK' num2str(number_mice(i))];
    end
end
genlist = [... % TEST =1, CTL =2
    1 1 1 1 ...
    1 2 2 2
    ];
genNames = {'Test', 'CTL'};
explist = [...% 50-50, 1
    2 2 2 2 ...
    2 2 2 2% 90-10, 2
   % 10-90, 3
   % 50-50, opto started after acquisition, 4
   % 90-10 VC, 1
   % 90-10, opto starts after learning, 2
    ];
expNames = {'90-10','90-10'}; % num of times to write exp names = num of genlist
struclist = [... % AC = 1, VC =2
    1 1 1 1 ...
    1 1 1 1];
strucNames = {'AC'};

plot_indiv_data = false;
plot_lick_psth = false;
savefig = true;
nFiles = 16; % = nb of days of behavior

% Initialization of variables
rdprime = nan(nSubj,nFiles);
pdprime = nan(nSubj,nFiles);
odprime = nan(nSubj,nFiles);
rfdprime = nan(nSubj,nFiles);
pdprime1 = nan(nSubj,nFiles);
pdprime2 = nan(nSubj,nFiles);
rcriterion = nan(nSubj,nFiles);
ocriterion = nan(nSubj,nFiles);
pcriterion = nan(nSubj,nFiles);
subjrates = cell(nSubj,1); 
subjnbctxtoutcome = cell(nSubj,1); 
mr_hit = nan(nSubj,nFiles); sr_hit = nan(nSubj,nFiles); % for lick latency
mp_hit = nan(nSubj,nFiles); sp_hit = nan(nSubj,nFiles);
mo_hit = nan(nSubj,nFiles); so_hit = nan(nSubj,nFiles);
mr_fa = nan(nSubj,nFiles); sr_fa = nan(nSubj,nFiles); % for lick latency
mp_fa = nan(nSubj,nFiles); sp_fa = nan(nSubj,nFiles);
mo_fa = nan(nSubj,nFiles); so_fa = nan(nSubj,nFiles);
mlr_hit = nan(nSubj,nFiles); slr_hit = nan(nSubj,nFiles); % for lick rate
mlp_hit = nan(nSubj,nFiles); slp_hit = nan(nSubj,nFiles);
mlo_hit = nan(nSubj,nFiles); slo_hit = nan(nSubj,nFiles);
mlr_fa = nan(nSubj,nFiles); slr_fa = nan(nSubj,nFiles); % for lick rate
mlp_fa = nan(nSubj,nFiles); slp_fa = nan(nSubj,nFiles);
mlo_fa = nan(nSubj,nFiles); slo_fa = nan(nSubj,nFiles);
lickhistcOFF = cell(nSubj,1); 
lickhistcON = cell(nSubj,1); 
lickhistrhit = cell(nSubj,1); 
lickhistrfa = cell(nSubj,1); 
lickhistohit = cell(nSubj,1); 
lickhistofa = cell(nSubj,1); 
lickhistphit = cell(nSubj,1); 
lickhistpfa = cell(nSubj,1); 
lickWindow = 1;
MAT = cell(nSubj,1);
%%
for nbsubj = 1:nSubj % through animals
    % Localize data
    subj = subjlist{nbsubj}; 
    if (str2double(subj(end-1:end))>=67 && str2double(subj(end-1:end))<=72) ||...
            str2double(subj(end-1:end))>76 || str2double(subj(end-2:end))>=100
        path = 'T:\su\DATA\behaviorData\opto\cohort1\'; % data 16
    else
        path = 'T:\su\DATA\behaviorData\opto\cohort1\'; % data 11
    end
    
%     if str2double(subj(end-2:end)) == 123
%         subj = '217';
%     elseif str2double(subj(end-2:end)) == 124
%         subj = '218';
%     elseif str2double(subj(end-2:end)) == 125
%         subj = '219';
%     elseif str2double(subj(end-2:end)) == 126
%         subj = '220';
%     elseif str2double(subj(end-2:end)) == 128
%         subj = '222';
%     end
%     if str2double(subj(end-2:end))>100
%         subj = lower(subj);
%     end  
    if nbsubj > 30 || explist(nbsubj)>= 10
        subjPath = [path subj '\GNG_LaserSineWavePhasicTone_Celine_ZZAW\Session Data\']; 
    else
        subjPath = [path subj '\GNG_LaserSineWavePhasicTone_Celine_ZZAW\Session Data\'];
    end
    files = dir([subjPath '*.mat']); 
    nFiles = length(files);
    if nFiles == 0
    end
    
    % Load and organize data
    filenames = cell(nFiles,1); 
    licks = cell(nFiles,1);
    matrix = []; 
    SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
    toremove = 0;
    for i=1:nFiles
        filenames{i} = [subjPath files(i).name];
        load(filenames{i}); 
        
        tempMat = nan(SessionData.nTrials,9);
        licks{i} = [];
        if i>1 && strcmp(files(i).date(1:11),files(i-1).date(1:11))  
            toremove = toremove+1;
        end
            
        numSess = i - toremove;
        tempMat(:,[SESS CTXT]) = [numSess*ones(SessionData.nTrials,1) SessionData.TrialSettings(1).context(1:SessionData.nTrials)]; 
        tempMat(:,START) = SessionData.TrialStartTimestamp;
        tempMat(:,STOP) = SessionData.TrialEndTimestamp;
        templicks = [];
        for j=1:SessionData.nTrials 
            tempMat(j,TONE_T) = SessionData.RawEvents.Trial{1,j}.States.Stimulus(1);
            if isfield(SessionData.RawEvents.Trial{1,j}.Events,'Port1In')
                l = SessionData.RawEvents.Trial{1,j}.Events.Port1In';
                templicks = [templicks;l + tempMat(j,START)];
                if ~isempty(find(l > tempMat(j,TONE_T)+0.1,1))
                    tempMat(j,LICKL) = l(find(l > tempMat(j,TONE_T)+0.1,1)) - (tempMat(j,TONE_T)+0.1);
                    tempMat(j,LICKR) = sum(l > tempMat(j,TONE_T)+0.1 &  l < tempMat(j,TONE_T)+0.1+lickWindow)/lickWindow;
                end
            elseif isfield(SessionData.RawEvents.Trial{1,j}.Events,'Port2In')
                l = SessionData.RawEvents.Trial{1,j}.Events.Port2In';
                templicks = [templicks;l + tempMat(j,START)];
                if ~isempty(find(l > tempMat(j,TONE_T)+0.1,1))
                    tempMat(j,LICKL) = l(find(l > tempMat(j,TONE_T)+0.1,1)) - (tempMat(j,TONE_T)+0.1);
                    tempMat(j,LICKR) = sum(l > tempMat(j,TONE_T)+0.1 &  l < tempMat(j,TONE_T)+0.1+lickWindow)/lickWindow;
                end
            end                
            if tempMat(j,CTXT) == 1 % opto
                if ~isnan(SessionData.RawEvents.Trial{1,j}.States.DrinkingLightON)
                    tempMat(j,OUTCOME) = 1; tempMat(j,TONE) = 1;
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.Miss)
                    tempMat(j,OUTCOME) = 2; tempMat(j,TONE) = 1; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.PunishLightON)
                    tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 2; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.CorrectReject) 
                    tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 2;
                end
            elseif tempMat(j,CTXT) == 2 || tempMat(j,CTXT) == 0 % reinforced, probe
                if ~isnan(SessionData.RawEvents.Trial{1,j}.States.Drinking)
                    tempMat(j,OUTCOME) = 1; tempMat(j,TONE) = 1;
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.Miss) 
                    tempMat(j,OUTCOME) = 2; tempMat(j,TONE) = 1; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.Punish)
                    tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 2; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.CorrectReject)
                    tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 2;
                end
            elseif tempMat(j,CTXT) == 4 % catch, light on
                if isfield(SessionData.RawEvents.Trial{1,j}.Events,'Port2In') % for bottom box, where licks were not detected
                    if any(sum(l > tempMat(j,TONE_T)+0.1+0.2 & l < tempMat(j,TONE_T)+0.1+0.2+2.5)) % if lick
                        tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 0;
                    else
                        tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 0; 
                    end
                else
                    if ~isnan(SessionData.RawEvents.Trial{1,j}.States.PunishLightON)
                        tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 0; 
                    elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.CorrectReject)
                        tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 0;
                    end
                end
            elseif tempMat(j,CTXT) == 3 % catch, light off
                if isfield(SessionData.RawEvents.Trial{1,j}.Events,'Port2In') % for bottom box, where licks were not detected
                    if any(sum(l > tempMat(j,TONE_T)+0.1+0.2 & l < tempMat(j,TONE_T)+0.1+0.2+2.5)) % if lick
                        tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 0;
                    else
                        tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 0; 
                    end
                else
                    if ~isnan(SessionData.RawEvents.Trial{1,j}.States.Punish)
                        tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 0; 
                    elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.CorrectReject)
                        tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 0;
                    end
                end
            end
        end 
        licks{numSess} = [licks{numSess};templicks];
        matrix = [matrix;tempMat];
    end
    
    ndays = max(matrix(:,SESS));
    rates = nan(ndays,14);
    nctxt = nan(ndays,14);
    nbins = 50;
    bins = linspace(-1,4,nbins);    
    lickhistr = nan(ndays,nbins);
    lickhisto = nan(ndays,nbins);
    lickhistcoff = nan(ndays,nbins);
    lickhistcon = nan(ndays,nbins);    
    lickhistr_hit = nan(ndays,nbins);    
    lickhistr_fa = nan(ndays,nbins);   
    lickhisto_hit = nan(ndays,nbins);   
    lickhisto_fa = nan(ndays,nbins);  
    lickhistp_hit = nan(ndays,nbins);   
    lickhistp_fa = nan(ndays,nbins); 
    for i=1:ndays
        templicks = licks{i};
        tempMat = matrix(matrix(:,SESS)==i,:);
        
        probeIdx = find(tempMat(:,CTXT)==0); 
        probeIdx1 = zeros(size(tempMat,1),1);
        probeIdx1(probeIdx(1):probeIdx(10)) = 1;
        probeIdx1 = logical(probeIdx1);
        probeIdx2 = zeros(size(tempMat,1),1);
        probeIdx2(probeIdx(11):probeIdx(end)) = 1;
        probeIdx2 = logical(probeIdx2);
    
        rhit = sum(tempMat(:,OUTCOME)==1 & tempMat(:,CTXT)==2) / sum(tempMat(:,CTXT)==2 & tempMat(:,TONE)==1); % hit in reinforced light off
        phit = sum(tempMat(:,OUTCOME)==1 & tempMat(:,CTXT)==0) / sum(tempMat(:,CTXT)==0 & tempMat(:,TONE)==1); % hit in probe
        ohit = sum(tempMat(:,OUTCOME)==1 & tempMat(:,CTXT)==1) / sum(tempMat(:,CTXT)==1 & tempMat(:,TONE)==1); % hit in opto (i.e. reinforced light on)
        rfhit = sum(tempMat(:,OUTCOME)==1 & (tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2)) / sum((tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2) & tempMat(:,TONE)==1); % hit in reinforced light off + on
        phit2_1 = sum(tempMat(:,OUTCOME)==1 & probeIdx1) / sum(probeIdx1 & tempMat(:,TONE)==1); % hit in probe, dividing probe bloc in 2
        phit2_2 = sum(tempMat(:,OUTCOME)==1 & probeIdx2) / sum(probeIdx2 & tempMat(:,TONE)==1); % hit in probe, dividing probe bloc in 2
        phit2 = [phit2_1 phit2_2];
        
        rfa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==2) / sum(tempMat(:,CTXT)==2 & tempMat(:,TONE)==2);
        pfa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==0) / sum(tempMat(:,CTXT)==0 & tempMat(:,TONE)==2);
        ofa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==1) / sum(tempMat(:,CTXT)==1 & tempMat(:,TONE)==2);
        rffa = sum(tempMat(:,OUTCOME)==3 & (tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2)) / sum((tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2) & tempMat(:,TONE)==2);
        pfa2_1 = sum(tempMat(:,OUTCOME)==3 & probeIdx1) / sum(probeIdx1 & tempMat(:,TONE)==2);
        pfa2_2 = sum(tempMat(:,OUTCOME)==3 & probeIdx2) / sum(probeIdx2 & tempMat(:,TONE)==2);
        pfa2 = [pfa2_1 pfa2_2];
        
        cofffa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==3) / sum(tempMat(:,CTXT)==3 & tempMat(:,TONE)==0); % catch light off FA
        conffa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==4) / sum(tempMat(:,CTXT)==4 & tempMat(:,TONE)==0); % catch light on FA
        rates(i,:) = [rhit rfa phit pfa ohit ofa rfhit rffa cofffa conffa phit2 pfa2];

        nctxt(i,:) = [sum(tempMat(:,CTXT)==2 & tempMat(:,TONE)==1) sum(tempMat(:,CTXT)==2 & tempMat(:,TONE)==2)... % numbers: e.g nb of HIT and nb of Target trials
            sum(tempMat(:,CTXT)==0 & tempMat(:,TONE)==1) sum(tempMat(:,CTXT)==0 & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)==1 & tempMat(:,TONE)==2) sum(tempMat(:,CTXT)==1 & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)~=0 & tempMat(:,TONE)==1) sum((tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2) & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)==3 & tempMat(:,TONE)==0) sum(tempMat(:,CTXT)==4 & tempMat(:,TONE)==0)...
            sum(probeIdx1 & tempMat(:,TONE)==1) sum(probeIdx2 & tempMat(:,TONE)==1)...
            sum(probeIdx2 & tempMat(:,TONE)==2) sum(probeIdx2 & tempMat(:,TONE)==2)];

        reinf = tempMat(:,CTXT)==2; opto = tempMat(:,CTXT)==1; probe = tempMat(:,CTXT)==0;
        coff = tempMat(:,CTXT)==3; con = tempMat(:,CTXT)==4;
        rt_r = RelativeTimes(templicks,[tempMat(reinf,START)+tempMat(reinf,TONE_T)-1 tempMat(reinf,START)+tempMat(reinf,TONE_T) tempMat(reinf,START)+tempMat(reinf,TONE_T)+4],[-1 0 4]);
        rt_o = RelativeTimes(templicks,[tempMat(opto,START)+tempMat(opto,TONE_T)-1 tempMat(opto,START)+tempMat(opto,TONE_T) tempMat(opto,START)+tempMat(opto,TONE_T)+4],[-1 0 4]);
        lickhistr(i,:) = hist(rt_r,bins)/sum(reinf);
        lickhisto(i,:) = hist(rt_o,bins)/sum(opto);
        hit = tempMat(:,OUTCOME)==1;fa = tempMat(:,OUTCOME)==3;
        rt_r_hit = RelativeTimes(templicks,[tempMat(reinf&hit,START)+tempMat(reinf&hit,TONE_T)-1 tempMat(reinf&hit,START)+tempMat(reinf&hit,TONE_T) tempMat(reinf&hit,START)+tempMat(reinf&hit,TONE_T)+4],[-1 0 4]);
        rt_o_hit = RelativeTimes(templicks,[tempMat(opto&hit,START)+tempMat(opto&hit,TONE_T)-1 tempMat(opto&hit,START)+tempMat(opto&hit,TONE_T) tempMat(opto&hit,START)+tempMat(opto&hit,TONE_T)+4],[-1 0 4]);
        rt_r_fa = RelativeTimes(templicks,[tempMat(reinf&fa,START)+tempMat(reinf&fa,TONE_T)-1 tempMat(reinf&fa,START)+tempMat(reinf&fa,TONE_T) tempMat(reinf&fa,START)+tempMat(reinf&fa,TONE_T)+4],[-1 0 4]);
        rt_o_fa = RelativeTimes(templicks,[tempMat(opto&fa,START)+tempMat(opto&fa,TONE_T)-1 tempMat(opto&fa,START)+tempMat(opto&fa,TONE_T) tempMat(opto&fa,START)+tempMat(opto&fa,TONE_T)+4],[-1 0 4]);
        lickhistr_hit(i,:) = hist(rt_r_hit,bins)/sum(reinf&hit);
        lickhistr_fa(i,:) = hist(rt_r_fa,bins)/sum(reinf&fa);
        lickhisto_hit(i,:) = hist(rt_o_hit,bins)/sum(opto&hit);
        lickhisto_fa(i,:) = hist(rt_o_fa,bins)/sum(opto&fa);
        rt_p_hit = RelativeTimes(templicks,[tempMat(probe&hit,START)+tempMat(probe&hit,TONE_T)-1 tempMat(probe&hit,START)+tempMat(probe&hit,TONE_T) tempMat(probe&hit,START)+tempMat(probe&hit,TONE_T)+4],[-1 0 4]);
        rt_p_fa = RelativeTimes(templicks,[tempMat(probe&fa,START)+tempMat(probe&fa,TONE_T)-1 tempMat(probe&fa,START)+tempMat(probe&fa,TONE_T) tempMat(probe&fa,START)+tempMat(probe&fa,TONE_T)+4],[-1 0 4]);
        lickhistp_hit(i,:) = hist(rt_p_hit,bins)/sum(probe&fa);
        lickhistp_fa(i,:) = hist(rt_p_fa,bins)/sum(probe&fa);
        
        rt_coff = RelativeTimes(templicks,[tempMat(coff,START)+tempMat(coff,TONE_T)-1 tempMat(coff,START)+tempMat(coff,TONE_T) tempMat(coff,START)+tempMat(coff,TONE_T)+4],[-1 0 4]);
        rt_con = RelativeTimes(templicks,[tempMat(con,START)+tempMat(con,TONE_T)-1 tempMat(con,START)+tempMat(con,TONE_T) tempMat(con,START)+tempMat(con,TONE_T)+4],[-1 0 4]);
        lickhistcoff(i,:) = hist(rt_coff,bins)/sum(coff);
        lickhistcon(i,:) = hist(rt_con,bins)/sum(con);
    end   
    
    lickhistcOFF{nbsubj} = lickhistcoff;
    lickhistcON{nbsubj} = lickhistcon;
    lickhistrhit{nbsubj} = lickhistr_hit;
    lickhistrfa{nbsubj} = lickhistr_fa;
    lickhistohit{nbsubj} = lickhisto_hit;
    lickhistofa{nbsubj} = lickhisto_fa;
    lickhistphit{nbsubj} = lickhistp_hit;
    lickhistpfa{nbsubj} = lickhistp_fa;
    
    ratescorr = rates;
    for j=1:size(rates,2)
        ratescorr(ratescorr(:,j)==0,j) = 1./(2*nctxt(ratescorr(:,j)==0,j)); 
        ratescorr(ratescorr(:,j)==1,j) = 1-1./(2*nctxt(ratescorr(:,j)==1,j)); 
    end
    
    subjrates{nbsubj} = ratescorr;
    subjnbctxtoutcome{nbsubj} = nctxt;
    
    rdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,1)) - norminv(ratescorr(:,2)) )'; % reinforced light off
    pdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,3)) - norminv(ratescorr(:,4)) )'; % probe
    odprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,5)) - norminv(ratescorr(:,6)) )'; % opto (i.e. reinforced light on)
    rfdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,7)) - norminv(ratescorr(:,8)) )'; % all reinforced (light on + off)
    rcriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,1)) + norminv(ratescorr(:,2)) )'/2; % reinforced light off
    ocriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,5)) + norminv(ratescorr(:,6)) )'/2; % opto (i.e. reinforced light on)
    pcriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,3)) + norminv(ratescorr(:,4)) )'/2; % probe
    
    pdprime1(nbsubj,1:ndays) = ( norminv(ratescorr(:,11)) - norminv(ratescorr(:,12)) )'; % probe 10 1st trials
    pdprime2(nbsubj,1:ndays) = ( norminv(ratescorr(:,13)) - norminv(ratescorr(:,14)) )'; % probe 10 last trials
    
    for i=1:nFiles
        mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL));
        sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL));
        mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        mo_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        so_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        
        mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        mo_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        so_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        
        mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        
        mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
    end
    MAT{nbsubj,1} = matrix;
    
    if plot_indiv_data
        % Plot subj data: H and FA rates, dprimes, lick latency, lick rate
        fig = figure; hold on;    

        subplot(2,4,1); hold on; % H and FA rates
        plot(1:ndays,rates(:,1),'k.-'); % hit reinf
        plot(1:ndays,rates(:,2),'k.--');% fa reinf
        plot(1:ndays,rates(:,3),'.-','color',[0.8 0.8 0.8]); % hit probe
        plot(1:ndays,rates(:,4),'.--','color',[0.8 0.8 0.8]); % fa probe
        plot(1:ndays,rates(:,5),'b.-'); % hit opto
        plot(1:ndays,rates(:,6),'b.--'); % fa opto
%         title([subjlist{nbsubj} ' (' strucNames{struclist(nbsubj)} ' ' expNames{explist(nbsubj)} ...
%             ', ' genNames{genlist(nbsubj)} ') - H and FA rates']);
%         legend('HR', 'FAR', 'HP', 'FAP', 'HO', 'FAO');

        title([subjlist{nbsubj}...
            ', ' genNames{genlist(nbsubj)} ' - H and FA rates']);
        ylim([0 1]);
        ylabel('Proportion');xlabel('Days');

        subplot(2,4,2); hold on; % d primes
        plot(1:ndays,rdprime(nbsubj,1:ndays),'k.-');
        plot(1:ndays,pdprime(nbsubj,1:ndays),'.-','color',[0.8 0.8 0.8]);
        plot(1:ndays,odprime(nbsubj,1:ndays),'b.-');    
        ylabel("d'");xlabel('Days');
        ylim([-1 5]);
%         legend('Reinforced','Probe','Opto');
        PlotHVLines(0,'h','k');
    
        subplot(2,4,5); hold on; % lick latency
        shadedErrorBar(1:ndays,mr_hit(nbsubj,1:ndays),sr_hit(nbsubj,1:ndays),{'color','k'},0.5);
        shadedErrorBar(1:ndays,mp_hit(nbsubj,1:ndays),sp_hit(nbsubj,1:ndays),{'color',[0.8 0.8 0.8]},0.5);
        shadedErrorBar(1:ndays,mo_hit(nbsubj,1:ndays),so_hit(nbsubj,1:ndays),{'color','b'},0.5);
        ylim([0 2.5]);ylabel('HIT Lick latency (1st lick) (s)');
        xlabel('Days');
        
        subplot(2,4,6); hold on; % lick latency
        shadedErrorBar(1:ndays,mr_fa(nbsubj,1:ndays),sr_fa(nbsubj,1:ndays),{'color','k'},0.5);
        shadedErrorBar(1:ndays,mp_fa(nbsubj,1:ndays),sp_fa(nbsubj,1:ndays),{'color',[0.8 0.8 0.8]},0.5);
        shadedErrorBar(1:ndays,mo_fa(nbsubj,1:ndays),so_fa(nbsubj,1:ndays),{'color','b'},0.5);
        ylim([0 2.5]);ylabel('FA Lick latency (1st lick) (s)');
        xlabel('Days');
        
        subplot(2,4,7); hold on; % lick rate FA
        shadedErrorBar(1:ndays,mlr_hit(nbsubj,1:ndays),slr_hit(nbsubj,1:ndays),{'color','k'},0.5);
        shadedErrorBar(1:ndays,mlp_hit(nbsubj,1:ndays),slp_hit(nbsubj,1:ndays),{'color',[0.8 0.8 0.8]},0.5);
        shadedErrorBar(1:ndays,mlo_hit(nbsubj,1:ndays),slo_hit(nbsubj,1:ndays),{'color','b'},0.5);
        ylim([0 10]);ylabel('HIT Lick rate (Hz)');
        xlabel('Days');
        
        subplot(2,4,8); hold on; % lick rate FA
        shadedErrorBar(1:ndays,mlr_fa(nbsubj,1:ndays),slr_fa(nbsubj,1:ndays),{'color','k'},0.5);
        shadedErrorBar(1:ndays,mlp_fa(nbsubj,1:ndays),slp_fa(nbsubj,1:ndays),{'color',[0.8 0.8 0.8]},0.5);
        shadedErrorBar(1:ndays,mlo_fa(nbsubj,1:ndays),slo_fa(nbsubj,1:ndays),{'color','b'},0.5);
        ylim([0 10]);ylabel('FA Lick rate (Hz)');
        xlabel('Days'); 
        fig.Position = [100 100 900 600];
        if savefig
            cd(pathsave);
%             saveas(fig,[subjlist{nbsubj} '-' ...
%                 strucNames{struclist(nbsubj)} '_' expNames{explist(nbsubj)} '_' ...
%                 genNames{genlist(nbsubj)} '-ActionAndPerformance.pdf']);
            saveas(fig,[subjlist{nbsubj} '-' genNames{genlist(nbsubj)} '-ActionAndPerformance.fig']);
            saveas(fig,[subjlist{nbsubj} '-' ...
                genNames{genlist(nbsubj)} '-ActionAndPerformance.pdf']);close(fig);
        end
        frow=ndays/7;
        fcol=ndays/frow;
        if plot_lick_psth
            fig=figure;hold on;
            for d=1:ndays
                subplot(frow,fcol,d);hold on;
                plot(bins,lickhistr_hit(d,:),'k');
                plot(bins,lickhistr_fa(d,:),'k:');
                plot(bins,lickhisto_hit(d,:),'b');
                plot(bins,lickhisto_fa(d,:),'b:');
                PlotHVLines(0,'v','k');
                title([subjlist{nbsubj} ', Day ' num2str(d)]);
                xlabel('Time from tone onset (s)');
            end
            fig.Position = [100 100 1200 500];
            if savefig
                cd(pathsave);
                saveas(fig,[subjlist{nbsubj} '-LickPSTHreinforced-'...
                    num2str(nbins) 'bins.fig']);
                saveas(fig,[subjlist{nbsubj} '-LickPSTHreinforced-' ...
                    num2str(nbins) 'bins.pdf']);
                close(fig);
            end

            fig=figure;hold on;
            for d=1:ndays
                subplot(frow,fcol,d);hold on;
                plot(bins,lickhistp_hit(d,:),'k');
                plot(bins,lickhistp_fa(d,:),'k:');
                PlotHVLines(0,'v','k');
                title([subjlist{nbsubj} ', Probe day ' num2str(d)]);
                xlabel('Time from tone onset (s)');
            end
            fig.Position = [100 100 1200 500];
            if savefig
                cd(pathsave);
                saveas(fig,[subjlist{nbsubj} '-LickPSTHprobe-' num2str(nbins) 'bins.fig']);
                saveas(fig,[subjlist{nbsubj} '-LickPSTHprobe-' ...
                    num2str(nbins) 'bins.pdf']);
                close(fig);
            end
        end
        
        drawnow;
    end        
end

%% Pool out variables to analyze

maxdays = nFiles;
rhit = nan(nSubj,maxdays);
rfa = nan(nSubj,maxdays);
phit = nan(nSubj,maxdays);
pfa = nan(nSubj,maxdays);
ohit = nan(nSubj,maxdays);
ofa = nan(nSubj,maxdays);
cofffa = nan(nSubj,maxdays);
confa = nan(nSubj,maxdays);

for s=1:nSubj
    if isempty(subjrates{s}), continue, end
    n = size(subjrates{s},1);
    
    rhit(s,1:n) = subjrates{s}(:,1);
    rfa(s,1:n) = subjrates{s}(:,2);
    
    phit(s,1:n) = subjrates{s}(:,3);
    pfa(s,1:n) = subjrates{s}(:,4);
    
    ohit(s,1:n) = subjrates{s}(:,5);
    ofa(s,1:n) = subjrates{s}(:,6);
    
    cofffa(s,1:n) = subjrates{s}(:,9);
    confa(s,1:n) = subjrates{s}(:,10);
end

%% Clean up variables (defined non stable / non learner)

status = ones(nSubj,1);
nonLearner = [];
acquiOnly = [];
[~,w_nonLearner] = ismember(nonLearner,number_mice);
status(w_nonLearner) = 0;
[~,w_acquiOnly] = ismember(acquiOnly,number_mice);
status(w_acquiOnly) = 2;
%%
unstable_subj = []; % mouse number
unstable_day =  []; % day until which the behavior is stable (remove from day+1)

for i=1:length(unstable_subj)
    subj = unstable_subj(i);
    fromday = unstable_day(i)+1;
    [~,w_subj] = ismember(subj,number_mice);
    rhit(w_subj,fromday:end) = nan;
    rfa(w_subj,fromday:end) = nan;
    phit(w_subj,fromday:end) = nan;
    pfa(w_subj,fromday:end) = nan;
    ohit(w_subj,fromday:end) = nan;
    ofa(w_subj,fromday:end) = nan;
    cofffa(w_subj,fromday:end) = nan;
    confa(w_subj,fromday:end) = nan;
    rdprime(w_subj,fromday:end) = nan;
    pdprime(w_subj,fromday:end) = nan;
    odprime(w_subj,fromday:end) = nan;
    rcriterion(w_subj,fromday:end) = nan;
    ocriterion(w_subj,fromday:end) = nan;
    pcriterion(w_subj,fromday:end) = nan;
    mr_hit(w_subj,fromday:end) = nan;
    mp_hit(w_subj,fromday:end) = nan;
    mo_hit(w_subj,fromday:end) = nan;
    mr_fa(w_subj,fromday:end) = nan;
    mp_fa(w_subj,fromday:end) = nan;
    mo_fa(w_subj,fromday:end) = nan;
    mlr_hit(w_subj,fromday:end) = nan;
    mlp_hit(w_subj,fromday:end) = nan;
    mlo_hit(w_subj,fromday:end) = nan;
    mlr_fa(w_subj,fromday:end) = nan;
    mlp_fa(w_subj,fromday:end) = nan;
    mlo_fa(w_subj,fromday:end) = nan;
end
%         
%% Compute percent correct
rpc = (rhit+(1-rfa))/2*100; 
opc = (ohit+(1-ofa))/2*100; 
ppc = (phit+(1-pfa))/2*100; 

%% Let's correct d' and percent correct in probe: fix after max
pdrimeOrignial = pdprime;
pdprime = pdrimeOrignial;
[~,wm] = max(pdprime(:,1:10),[],2);
for i=1:nSubj
    pdprime(i,wm(i):end) = pdprime(i,wm(i));
end    

ppcOrignial = ppc;
[~,wm] = max(ppc(:,1:10),[],2);
for i=1:nSubj
    ppc(i,wm(i):end) = ppc(i,wm(i));
end    
%% Analysis by experiment

VAR = {rhit ohit phit; rfa ofa pfa; rdprime odprime pdprime; rpc opc ppc;...
    rcriterion ocriterion pcriterion};
ylims = [[0 1];[0 1];[-1 5];[40 100];[-2 2]];
nExp = numel(unique(explist));
condNames = {'Reinf OFF','Reinf ON','Probe'};
varNames = {'hit','fa','dprime','PercentCorrect','Criterion'};
savefig = true;

% for i=1:nExp
%     ok = explist==i;
%     nStruc = numel(unique(struclist(ok)));
%     for j=1:nStruc
%         ok = explist==i & struclist==j;
%         for v=1:size(var,2)             
%             fig = figure;
%             set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
%             var = VAR(v,:);
%             for c=1:3 % loop trough conditions, i.e. reinf, opto, probe
%                 if c==3 
%                     taken = (status==1 | status==2)';
%                 else
%                     taken = (status==1)';
%                 end
%                 groups = genlist(ok & taken);            
%                 data = var{c}(ok & taken,:);
%                 subplot(2,size(var,2),c);hold on;
%                 plot(data(groups==2,:)','k');
%                 plot(data(groups==1,:)','c');
%                 title([expNames{i} '-' strucNames{j} ', ' condNames{c} ' (ctl=' num2str(sum(groups==2))...
%                     ', test=' num2str(sum(groups==1)) ')']);
%                 ylabel(varNames{v});
%                 xlabel('Days');
%                 ylim(ylims(v,:));
%             end
%         end        
%     end
% end

for i=1:nExp
    ok = explist==2;
    nStruc = numel(unique(struclist(ok)));
    for j=1:nStruc
%         ok = explist==i & struclist==j;
        ok=explist==2;
        for v=1:size(VAR,1)             
            fig = figure;
            set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
            orient(fig,'landscape')
            var = VAR(v,:);
            for c=1:3 % loop trough conditions, i.e. reinf, opto, probe
                if c==3 
                    taken = (status==1 | status==2)';
                else
                    taken = (status==1)';
                end
                groups = genlist(ok & taken);            
                data = var{c}(ok & taken,:);
                    data=var{c};
                subplot(2,size(var,2),c);hold on;
                plot(nanmean(data(groups==2,:))','k','linewidth',3);
                plot(nanmean(data(groups==2,:))'+nansem(data(groups==2,:))','k');
                plot(nanmean(data(groups==2,:))'-nansem(data(groups==2,:))','k');
                plot(nanmean(data(groups==1,:))','c','linewidth',3);
                plot(nanmean(data(groups==1,:))'+nansem(data(groups==1,:))','c');
                plot(nanmean(data(groups==1,:))'-nansem(data(groups==1,:))','c');
                title([expNames{i} '-' strucNames{j} ', ' ...
                    condNames{c} ' (ctl=' num2str(sum(groups==2))...                    
                    ', test=' num2str(sum(groups==1)) ')']);
                ylabel(varNames{v});
                xlabel('Days');
                ylim(ylims(v,:));
                xlim([0 30]);
                if v==5, PlotHVLines(0,'h','k'); end
            end
            if savefig
                cd(pathsave);
                t = [expNames{i} '-' strucNames{j} '-' varNames{v}];
                saveas(fig,[t '.pdf']);
                close(fig);
            end
        end        
    end
end

%% Stats: difference percent correct
v = 4;
var = VAR(v,:);
maxday = nFiles;
for i=1:nExp
    ok = explist==i;
    nStruc = numel(unique(struclist(ok)));
    for j=1:nStruc
        ok = explist==i & struclist==j;
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape')
        for c=1:3 % loop trough conditions, i.e. reinf, opto, probe
            if c==3 
                taken = (status==1 | status==2)';
            else
                taken = (status==1)';
            end
            groups = genlist(ok & taken);            
            data = var{c}(ok & taken,:);
            subplot(2,size(var,2),c);hold on;
            plot(nanmean(data(groups==2,1:maxday))','k','linewidth',3);
            plot(nanmean(data(groups==2,1:maxday))'+nansem(data(groups==2,1:maxday))','k');
            plot(nanmean(data(groups==2,1:maxday))'-nansem(data(groups==2,1:maxday))','k');
            plot(nanmean(data(groups==1,1:maxday))','c','linewidth',3);
            plot(nanmean(data(groups==1,1:maxday))'+nansem(data(groups==1,1:maxday))','c');
            plot(nanmean(data(groups==1,1:maxday))'-nansem(data(groups==1,1:maxday))','c');
            title([expNames{i} '-' strucNames{j} ', ' condNames{c} ' (ctl=' num2str(sum(groups==2))...
                ', test=' num2str(sum(groups==1)) ')']);
            ylabel(varNames{v});
            xlabel('Days');
            ylim(ylims(v,:));
            xlim([0 21]);
        end
    end        
end

%% Plot catch trials
savefig = false;
for i=1:nExp
%     ok = explist==i;
    ok=explist==i+1;
    nStruc = numel(unique(struclist(ok)));
    for j=1:nStruc
%         ok = explist==i & struclist==j;
        ok=explist==i+1;
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape');
        taken = (status==1)';        
        groups = genlist(ok & taken);   
        dataOFF = cofffa(ok & taken,:);   
        dataON = confa(ok & taken,:);  
        
        subplot(2,2,1);hold on;
        doff = dataOFF(groups==2,:)';
        don = dataON(groups==2,:)';            
        if ismember(i,[1:nExp])
            doff(22:end,:) = nan;don(22:end,:) = nan;
        end
        bar([1 2],nanmean([doff(:) don(:)])');
        plot([doff(:) don(:)]','k.-');
        ylabel('FA catch');ylim([0 1]);
        [~,p] = ttest(doff(:),don(:));
        set(gca,'xtick',[1 2],'xticklabel',{'OFF','ON'});
        title([expNames{i} '-' strucNames{j} ', CTL (n=' num2str(sum(groups==2)) ', ttestpaired=' num2str(p) ')']);

        subplot(2,2,2);hold on;
        doff = dataOFF(groups==1,:)';
        don = dataON(groups==1,:)';    
        if ismember(i,[1:nExp])
            doff(22,:) = nan;don(22,:) = nan;
        end
        bar([1 2],nanmean([doff(:) don(:)])');
        plot([doff(:) don(:)]','k.-');
        ylabel('FA catch');ylim([0 1]);
        [~,p] = ttest(doff(:),don(:));
        set(gca,'xtick',[1 2],'xticklabel',{'OFF','ON'});
        title([expNames{i} '-' strucNames{j} ', TEST (n=' num2str(sum(groups==2)) ', ttestpaired=' num2str(p) ')']);      
        
        if savefig
            cd(pathsave);
            t = [expNames{i} '-' strucNames{j} '-Catch_Trials'];
            saveas(fig,[t '.pdf']);
            close(fig);
        end
    end
end

%% Lick rate and lick latency

   VAR = {mr_hit mo_hit mp_hit; mr_fa mo_fa mp_fa; mlr_hit mlo_hit mlp_hit; mlr_fa mlo_fa mlp_fa};
   ylims = [[0 2.5];[0 2.5];[0 10];[0 10]];
  nExp = numel(unique(explist));
 condNames = {'Reinf OFF','Reinf ON','Probe'};
   varNames = {'HIT lick latency','FA lick latency','HIT lick rate','FA lick rate'};
   savefig = true;

   for i=1:nExp
%        ok = explist==i;
        ok=explist==i+1;
       nStruc = numel(unique(struclist(ok)));
       for j=1:nStruc
%            ok = explist==i & struclist==j;
            ok=explist==i+1;
           for v=1:size(VAR,1)             
               fig = figure;
               set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
               orient(fig,'landscape')
               var = VAR(v,:);
               for c=1:3 % loop trough conditions, i.e. reinf, opto, probe
                   if c==3 
                       taken = (status==1 | status==2)';
                   else
                       taken = (status==1)';
                   end
                   groups = genlist(ok & taken);            
                   data = var{c}(ok & taken,:);
                   subplot(2,size(var,2),c);hold on;
                   plot(nanmean(data(groups==2,:))','k','linewidth',3);
                   plot(nanmean(data(groups==2,:))'+nansem(data(groups==2,:))','k');
                   plot(nanmean(data(groups==2,:))'-nansem(data(groups==2,:))','k');
                   plot(nanmean(data(groups==1,:))','c','linewidth',3);
                   plot(nanmean(data(groups==1,:))'+nansem(data(groups==1,:))','c');
                   plot(nanmean(data(groups==1,:))'-nansem(data(groups==1,:))','c');
                   title([expNames{i} '-' strucNames{j} ', ' condNames{c} ' (ctl=' num2str(sum(groups==2))...
                       ', test=' num2str(sum(groups==1)) ')']);
                   ylabel(varNames{v});
                   xlabel('Days');
                   ylim(ylims(v,:));
                   xlim([0 30]);
                   if v==5, PlotHVLines(0,'h','k'); end
               end
               if savefig
                   cd(pathsave);
                   t = [expNames{i} '-' strucNames{j} '-' varNames{v}];
                   saveas(fig,[t '.pdf']);
                   close(fig);
               end
           end        
       end
   end

%% Plot difference performance light-off light-on

maxday = 16;
days = 1:maxday;
savefig = false; 
for i=1:nExp
%     ok = explist==i;
    ok = explist==i+1;
    nStruc = numel(unique(struclist(ok)));
    for j=1 % :nStruc
%         ok = explist==i & struclist==j;
        ok = explist==i+1;
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape');
        taken = (status==1)';        
        groups = genlist(ok & taken);            
        
        subplot(1,2,1);hold on;
        data = opc(ok & taken,days)-rpc(ok & taken,days);  
        plot(nanmean(data(groups==2,:))','k','linewidth',3);
        plot(nanmean(data(groups==2,:))'+nansem(data(groups==2,:))','k');
        plot(nanmean(data(groups==2,:))'-nansem(data(groups==2,:))','k');
        plot(nanmean(data(groups==1,:))','c','linewidth',3);
        plot(nanmean(data(groups==1,:))'+nansem(data(groups==1,:))','c');
        plot(nanmean(data(groups==1,:))'-nansem(data(groups==1,:))','c');
        PlotHVLines(0,'h','k');
        ylim([-20 10]);
        xlim([1 16]);
        ylabel('Diff percent correct');
        xlabel([0 22]);xlabel('Days');
        title([expNames{i} '-' strucNames{j} ' (ctl=' num2str(sum(groups==2)) ...
            ', test=' num2str(sum(groups==1)) ')']);      
        
        subplot(1,2,2);hold on;
        dataR = rpc(ok & taken,days);  
        dataO = opc(ok & taken,days);
        plot(dataR(groups==1,days),dataO(groups==1,days),'c.','markersize',10);
        plot(dataR(groups==2,days),dataO(groups==2,days),'k.','markersize',10);
        xlabel('Percent correct light-off');
        ylabel('Percent correct light-on');
        PlotHVLines(50,'v','k');
        PlotHVLines(50,'h','k');
        ylim([35 100]);xlim([35 100]);
        plot(35:100,35:100,'k');
        
        rdprime_test = dataR(groups==1,days);rdprime_test = rdprime_test(:);
        odprime_test = dataO(groups==1,days);odprime_test = odprime_test(:);
        nop = any(isnan([rdprime_test odprime_test]),2);
        [rho_test,yfit] = LinearRegression(rdprime_test(~nop),odprime_test(~nop));
        plot(rdprime_test(~nop),yfit,'c','linewidth',2);
        [~,p_test] = corr(rdprime_test(~nop),odprime_test(~nop)); 
        rdprime_ctl = dataR(groups==2,days);rdprime_ctl = rdprime_ctl(:);
        odprime_ctl = dataO(groups==2,days);odprime_ctl = odprime_ctl(:);
        nop = any(isnan([rdprime_ctl odprime_ctl]),2);
        [rho_ctl,yfit] = LinearRegression(rdprime_ctl(~nop),odprime_ctl(~nop)); 
        plot(rdprime_ctl(~nop),yfit,'k','linewidth',2);
        [~,p_ctl] = corr(rdprime_ctl(~nop),odprime_ctl(~nop)); 
        title(['r2 ctl = ' num2str(rho_ctl) ', test = ' num2str(rho_test)]);
        
        if savefig
            cd(pathsave);
            t = [expNames{i} '-' strucNames{j} '-Diff_Perf_ON_OFF'];
            saveas(fig,[t '.pdf']);
            close(fig);
        end
        
        if i==1 & j==1
            fig = figure;
            set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
            orient(fig,'landscape');
            subplot(2,2,1);hold on;
            dataR = rpc(ok & taken,days);  
            dataO = opc(ok & taken,days);
            plot(dataR(groups==2,days),dataO(groups==2,days),'k.','markersize',10);
%             totest = find(groups==1);
            lightdriven = zeros(1,length(groups));
%             lightdriven([7;8;9;10]) = 1;
            lightdriven([7;9;10;21]) = 1; % Need to manually change
            lightdriven = logical(lightdriven);
            plot(dataR(groups==1 & lightdriven,days),dataO(groups==1 & lightdriven,days)...
                ,'b.','markersize',10);         
            plot(dataR(groups==1 & ~lightdriven,days),dataO(groups==1 & ~lightdriven,days)...
                ,'c.','markersize',10);                 
            xlabel('Percent correct light-off');
            ylabel('Percent correct light-on');
            PlotHVLines(50,'v','k');
            PlotHVLines(50,'h','k');
            ylim([35 100]);xlim([35 100]);
            plot(35:100,35:100,'k');
            
            subplot(2,2,2);hold on;
            data = opc(ok & taken,days)-rpc(ok & taken,days);  
            plot(nanmean(data(groups==2,:))','k','linewidth',3);
            plot(nanmean(data(groups==2,:))'+nansem(data(groups==2,:))','k');
            plot(nanmean(data(groups==2,:))'-nansem(data(groups==2,:))','k');
            plot(nanmean(data(groups==1 & lightdriven,:))','b','linewidth',3);
            plot(nanmean(data(groups==1 & lightdriven,:))'+nansem(data(groups==1 & lightdriven,:))','b');
            plot(nanmean(data(groups==1 & lightdriven,:))'-nansem(data(groups==1 & lightdriven,:))','b');
            plot(nanmean(data(groups==1 & ~lightdriven,:))','c','linewidth',3);
            plot(nanmean(data(groups==1 & ~lightdriven,:))'+nansem(data(groups==1 & ~lightdriven,:))','c');
            plot(nanmean(data(groups==1 & ~lightdriven,:))'-nansem(data(groups==1 & ~lightdriven,:))','c');
            PlotHVLines(0,'h','k');
            ylim([-50 10]);
            ylabel('Diff percent correct');
            xlabel([0 22]);xlabel('Days');
            title([expNames{i} '-' strucNames{j} ' (ctl=' num2str(sum(groups==2)) ...
                ', lightdrivent=' num2str(sum(groups==1 & lightdriven)) ...
                ', not lightdrivent=' num2str(sum(groups==1 & ~lightdriven)) ')']);
            
            % plot catch light on vs catch light off
            dataOFF = cofffa(ok & taken,:);   
            dataON = confa(ok & taken,:);  
            
            subplot(2,2,3);hold on;        
            doff = dataOFF(groups==1 & lightdriven,:)';
            don = dataON(groups==1 & lightdriven,:)';            
            doff(22:end,:) = nan;don(22:end,:) = nan;
            bar([1 2],nanmean([doff(:) don(:)])');
            plot([doff(:) don(:)]','k.-');
            ylabel('FA catch');ylim([0 1]);
            [~,p] = ttest(doff(:),don(:));
            set(gca,'xtick',[1 2],'xticklabel',{'OFF','ON'});
            title([expNames{i} '-' strucNames{j} ...
                ', TEST lightdriven (n=' num2str(sum(groups==1 & lightdriven))...
                '), ttestpaired=' num2str(p) ')']);
            
            subplot(2,2,4);hold on;        
            doff = dataOFF(groups==1 & ~lightdriven,:)';
            don = dataON(groups==1 & ~lightdriven,:)';            
            doff(22:end,:) = nan;don(22:end,:) = nan;
            bar([1 2],nanmean([doff(:) don(:)])');
            plot([doff(:) don(:)]','k.-');
            ylabel('FA catch');ylim([0 1]);
            [~,p] = ttest(doff(:),don(:));
            set(gca,'xtick',[1 2],'xticklabel',{'OFF','ON'});
            title([expNames{i} '-' strucNames{j} ...
                ', TEST NOT lightdriven (n=' num2str(sum(groups==1 & ~lightdriven))...
                '), ttestpaired=' num2str(p) ')']);
            if savefig
                cd(pathsave);
                t = [expNames{i} '-' strucNames{j} '-Diff_Perf_ON_OFF_GroupLightDriven'];
                saveas(fig,[t '.pdf']);
                close(fig);
            end           
            
        elseif i==4
            fig = figure;
            set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
            orient(fig,'landscape');
            subplot(2,2,1);hold on;
            dataR = rpc(ok & taken,days);  
            dataO = opc(ok & taken,days);
            plot(dataR(groups==2,days),dataO(groups==2,days),'k.','markersize',10);
            lightdriven = zeros(1,length(groups));
            lightdriven([2;3;7]) = 1;
            lightdriven = logical(lightdriven);
            plot(dataR(groups==1 & lightdriven,days),dataO(groups==1 & lightdriven,days)...
                ,'b.','markersize',10);         
            plot(dataR(groups==1 & ~lightdriven,days),dataO(groups==1 & ~lightdriven,days)...
                ,'c.','markersize',10);                 
            xlabel('Percent correct light-off');
            ylabel('Percent correct light-on');
            PlotHVLines(50,'v','k');
            PlotHVLines(50,'h','k');
            ylim([35 100]);xlim([35 100]);
            plot(35:100,35:100,'k');
            
            subplot(2,2,2);hold on;
            data = opc(ok & taken,days)-rpc(ok & taken,days);  
            plot(nanmean(data(groups==2,:))','k','linewidth',3);
            plot(nanmean(data(groups==2,:))'+nansem(data(groups==2,:))','k');
            plot(nanmean(data(groups==2,:))'-nansem(data(groups==2,:))','k');
            plot(nanmean(data(groups==1 & lightdriven,:))','b','linewidth',3);
            plot(nanmean(data(groups==1 & lightdriven,:))'+nansem(data(groups==1 & lightdriven,:))','b');
            plot(nanmean(data(groups==1 & lightdriven,:))'-nansem(data(groups==1 & lightdriven,:))','b');
            plot(data(groups==1 & ~lightdriven,:)','c','linewidth',3);
            PlotHVLines(0,'h','k');
            ylim([-50 10]);
            ylabel('Diff percent correct');
            xlabel([0 22]);xlabel('Days');
            title([expNames{i} '-' strucNames{j} ' (ctl=' num2str(sum(groups==2)) ...
                ', lightdrivent=' num2str(sum(groups==1 & lightdriven)) ...
                ', not lightdrivent=' num2str(sum(groups==1 & ~lightdriven)) ')']);
            
            % plot catch light on vs catch light off
            dataOFF = cofffa(ok & taken,:);   
            dataON = confa(ok & taken,:);              
            subplot(2,2,3);hold on;cla            
            doff = dataOFF(groups==1 & lightdriven,:)';
            don = dataON(groups==1 & lightdriven,:)';            
            plot(doff,'k');
            plot(don,'b');
            xlabel('Days');
            ylabel('FA catch');ylim([0 1]);
            title([expNames{i} '-' strucNames{j} ...
                ', TEST lightdriven (n=' num2str(sum(groups==1 & lightdriven))...
                ')']);
            
            subplot(2,2,4);hold on;        
            doff = dataOFF(groups==1 & ~lightdriven,:)';
            don = dataON(groups==1 & ~lightdriven,:)';            
            plot(doff,'k');
            plot(don,'b');
            xlabel('Days');
            ylabel('FA catch');ylim([0 1]);
            set(gca,'xtick',[1 2],'xticklabel',{'OFF','ON'});
            title([expNames{i} '-' strucNames{j} ...
                ', TEST NOT lightdriven (n=' num2str(sum(groups==1 & ~lightdriven))...
                ')']);
            
            if savefig
                cd(pathsave);
                t = [expNames{i} '-' strucNames{j} '-Diff_Perf_ON_OFF_GroupLightDriven'];
                saveas(fig,[t '.pdf']);
                close(fig);
            end   
        end
    end
end

%% Lick PSTH
savefig = true;
for i=[1 4]
    ok = explist==i;
    nStruc = numel(unique(struclist(ok)));
    for j=1%:nStruc
        ok = explist==i & struclist==j;
        taken = (status==1)';
        groups = genlist(ok & taken);            
        datarhit = lickhistrhit(ok & taken,:);
        dataohit = lickhistohit(ok & taken,:);

        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape')
        drctl = datarhit(groups==2,:);
        doctl = dataohit(groups==2,:);
        ndays = min([size(drctl{1,1},1);size(drctl{end,1},1)]);
        nctl = size(drctl,1);
        mrlickdays = nan(ndays,50);
        molickdays = nan(ndays,50);
        for dd = 1:ndays 
            rlickdays = nan(nctl,50);
            olickdays = nan(nctl,50);
            for ct = 1:nctl
                rlickdays(ct,:) = drctl{ct}(dd,:);
                olickdays(ct,:) = doctl{ct}(dd,:);
            end
            mrlickdays(dd,:) = nanmean(rlickdays);
            molickdays(dd,:) = nanmean(olickdays);
            subplot(5,5,dd);hold on;
            plot(mrlickdays(dd,:),'k');hold on;
            plot(molickdays(dd,:),'b');
%             xlim([5 15]);
            title([expNames{i} '-D' num2str(dd) ', ctl']);
        end
        
        if savefig
            cd(pathsave);
            t = [expNames{i} '-' strucNames{j} '-LickPSTH-Ctl'];
            saveas(fig,[t '.pdf']);
            close(fig);
        end   
%         lightdriven = zeros(1,length(groups));
%         if i==1        
%             lightdriven([7;9;10]) = 1;
%         elseif i==4
%             lightdriven([2;3;7]) = 1;
%         end
%         lightdriven = logical(lightdriven);
        lightdriven = zeros(1,length(groups));
        if i==1        
            lightdriven([7;9;10]) = 1;
        elseif i==4
            lightdriven([2;6]) = 1;
        end
        lightdriven = logical(lightdriven);
            
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape')
        drtest = datarhit(groups==1 & lightdriven,:);
        dotest = dataohit(groups==1 & lightdriven,:);
        ndays = min([size(drtest{1,1},1);size(drtest{end,1},1)]);
        nctl = size(drtest,1);
        mrlickdays = nan(ndays,50);
        molickdays = nan(ndays,50);
        for dd = 1:ndays 
            rlickdays = nan(nctl,50);
            olickdays = nan(nctl,50);
            for ct = 1:nctl
                rlickdays(ct,:) = drtest{ct}(dd,:);
                olickdays(ct,:) = dotest{ct}(dd,:);
            end
            mrlickdays(dd,:) = nanmean(rlickdays);
            molickdays(dd,:) = nanmean(olickdays);
            subplot(5,5,dd);hold on;
            plot(mrlickdays(dd,:),'k');hold on;
            plot(molickdays(dd,:),'b');
%             xlim([5 15]);
            title([expNames{i} '-D' num2str(dd) ', light driven']);
        end
        
         if savefig
            cd(pathsave);
            t = [expNames{i} '-' strucNames{j} '-LickPSTH-LightDriven'];
            saveas(fig,[t '.pdf']);
            close(fig);
        end   
        
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape')
        drtest = datarhit(groups==1 & ~lightdriven,:);
        dotest = dataohit(groups==1 & ~lightdriven,:);
        ndays = size(drtest{end,1},1);
        nctl = size(drtest,1);
        mrlickdays = nan(ndays,50);
        molickdays = nan(ndays,50);
        for dd = 1:ndays 
            rlickdays = nan(nctl,50);
            olickdays = nan(nctl,50);
            for ct = 1:nctl
                rlickdays(ct,:) = drtest{ct}(dd,:);
                olickdays(ct,:) = dotest{ct}(dd,:);
            end
            mrlickdays(dd,:) = nanmean(rlickdays,1);
            molickdays(dd,:) = nanmean(olickdays,1);
            subplot(5,5,dd);hold on;
            plot(mrlickdays(dd,:),'k');hold on;
            plot(molickdays(dd,:),'b');
%             xlim([5 15]);
            title([expNames{i} '-D' num2str(dd) ', no light driven']);
        end
        
        if savefig
            cd(pathsave);
            t = [expNames{i} '-' strucNames{j} '-LickPSTH-NoLightDriven'];
            saveas(fig,[t '.pdf']);
            close(fig);
        end
    end
end
       
%% Lick rate and lick latency diff between ctl, light driven and no light driven

VAR = {mr_hit mo_hit mp_hit; mr_fa mo_fa mp_fa; mlr_hit mlo_hit mlp_hit; mlr_fa mlo_fa mlp_fa};
% ylims = [[-0.1 0.1];[-1.5 0.5];[-1 1];[-1 6]]; % diff
ylims = [[0 2];[0 2];[0 9];[0 9]];
varNames = {'HIT lick latency','FA lick latency','HIT lick rate','FA lick rate'};
savefig = false;

for i=[1 4]
    ok = explist==i;
    for j=1 %:nStruc
        ok = explist==i & struclist==j;
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape');
        for v=1:size(VAR,1)             
            
            var = VAR(v,:);            
            taken = (status==1)';

            groups = genlist(ok & taken);            
%             diffdata = var{2}(ok & taken,:)-var{1}(ok & taken,:);
            diffdata = var{2}(ok & taken,:);
            
            subplot(2,2,v);hold on;
            
            plot(nanmean(diffdata(groups==2,:))','k','linewidth',3);
            plot(nanmean(diffdata(groups==2,:))'+nansem(diffdata(groups==2,:))','k');
            plot(nanmean(diffdata(groups==2,:))'-nansem(diffdata(groups==2,:))','k');
            
            lightdriven = zeros(1,length(groups));
            if i==1        
                lightdriven([7;9;10]) = 1;
            elseif i==4
                lightdriven([2;6]) = 1;
            end
            lightdriven = logical(lightdriven);

            plot(nanmean(diffdata(groups==1 & lightdriven,:))','b','linewidth',3);
            plot(nanmean(diffdata(groups==1 & lightdriven,:))'+nansem(diffdata(groups==1 & lightdriven,:))','b');
            plot(nanmean(diffdata(groups==1 & lightdriven,:))'-nansem(diffdata(groups==1 & lightdriven,:))','b');
            
            
            plot(nanmean(diffdata(groups==1 & ~lightdriven,:),1)','c','linewidth',3);
            if sum(groups==1 & ~lightdriven)>1
                plot(nanmean(diffdata(groups==1 & ~lightdriven,:))'+nansem(diffdata(groups==1 & ~lightdriven,:))','c');
                plot(nanmean(diffdata(groups==1 & ~lightdriven,:))'-nansem(diffdata(groups==1 & ~lightdriven,:))','c');
            end
            title([expNames{i} '-' strucNames{j} ' (ctl=' num2str(sum(groups==2))...
                ', light=' num2str(sum(groups==1 & lightdriven)) ', '...
                ', notlight=' num2str(sum(groups==1 & ~lightdriven)) ')']);
                        
            ylabel(varNames{v});
            xlabel('Days');
            ylim(ylims(v,:));
            xlim([0 30]);
            PlotHVLines(0,'h','k'); 
        end        
        if savefig
            cd(pathsave);
            t = [expNames{i} '-' strucNames{j} '-DifferenceLickLightON-OFF'];
            saveas(fig,[t '.pdf']);
            close(fig);
        end
    end
end

%% Correlation difference catch LIGHT ON vs OFF and HIT lick latency 

for i=[1 4]
    ok = explist==i;
    for j=1 % only AC
        ok = explist==i & struclist==j;
        taken = (status==1)';
        groups = genlist(ok & taken);          
        latO = mo_fa(ok & taken,:);
        latR =  mr_fa(ok & taken,:);
        diff_latency = latR-latO; 
        dataOFF = cofffa(ok & taken,:);   
        dataON = confa(ok & taken,:); 
        diff_facatch = dataOFF-dataON;    
                
        % control        
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape');
        nctls = sum(groups==2);
        ctls = find(groups==2);
        for c = 1:nctls
            subplot(nctls,2,c+(c-1));hold on;
%             plot(diff_latency(ctls(c),:)','k');
            plot(latO(ctls(c),:)','b');
            xlim([0 20]);ylim([0 3]);
            PlotHVLines(0,'h','k--');
            subplot(nctls,2,c+1+(c-1));hold on;
            plot(dataOFF(ctls(c),:)','k');
            plot(dataON(ctls(c),:)','b');
            xlim([0 18]);
        end
                
        % tests        
        fig = figure;
        set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
        orient(fig,'landscape');
        ntests = sum(groups==1);
        tests = find(groups==1);
        for c = 1:ntests
            subplot(ntests,2,c+(c-1));hold on;
%             plot(diff_latency(tests(c),:)','k');
            plot(latO(tests(c),:)','b');
            xlim([0 20]);ylim([0 1]);
            PlotHVLines(0,'h','k--');
            subplot(ntests,2,c+1+(c-1));hold on;
            plot(dataOFF(tests(c),:)','k');
            plot(dataON(tests(c),:)','b');
            xlim([0 18]);
        end        
        
        % Scatter plot: FA latency Reinf light ON = f(FA catch Light ON);
        fig = figure;
        subplot(2,2,1);hold on;
        for c = 1:nctls
            d = find(~isnan(latO(ctls(c),:)'),1);
            maxd = size(dataON,2);
            plot(dataON(ctls(c),d:maxd)',latR(ctls(c),d:maxd)','k.');
        end
        ylim([0 3]);xlim([0 1]);
        ylabel('FA lick latency Reinf');
        xlabel('FA catch light ON');
        legend('Ctl');
        x = dataON(ctls,d:maxd)';x = x(:);
        y = latR(ctls,d:maxd)';y = y(:);
        nop = isnan(x) | isnan(y);
        [rho, yfit] = LinearRegression(x(~nop),y(~nop));
        plot(x(~nop),yfit,'k');
        
        subplot(2,2,2);hold on;
        for c = 1:ntests
            d = find(~isnan(latO(tests(c),:)'),1);
            maxd = size(dataON,2);
            plot(dataON(tests(c),d:maxd)',latR(tests(c),d:maxd)','c.');
        end
        ylim([0 3]);
        ylabel('FA lick latency Reinf');
        xlabel('FA catch light ON');
        legend('Test');
        
        lightdriven = zeros(1,length(groups));
        if i==1        
            lightdriven([7;9;10]) = 1;
        elseif i==4
            lightdriven([2;6]) = 1;
        end
        lightdriven = logical(lightdriven);
        
        subplot(2,2,3);hold on;
        ntestsD = sum(groups==1 & lightdriven);
        testsD = find(groups==1 & lightdriven);
        for c = 1:ntestsD
            d = find(~isnan(latO(testsD(c),:)'),1);
            maxd = size(dataON,2);
            plot(dataON(testsD(c),d:maxd)',latR(testsD(c),d:maxd)','b.');
        end
        x = dataON(testsD,d:maxd)';x = x(:);
        y = latR(testsD,d:maxd)';y = y(:);
        nop = isnan(x) | isnan(y);
        [rho, yfit] = LinearRegression(x(~nop),y(~nop));
        plot(x(~nop),yfit,'b');
        [c,p]= corr(x(~nop),y(~nop));
        
        ntestsND = sum(groups==1 & ~lightdriven);
        testsND = find(groups==1 & ~lightdriven);
        for c = 1:ntestsND
            d = find(~isnan(latO(testsND(c),:)'),1);
            maxd = size(dataON,2);
            plot(dataON(testsND(c),d:maxd)',latR(testsND(c),d:maxd)','c.');
        end
        x = dataON(testsND,d:maxd)';x = x(:);
        y = latR(testsND,d:maxd)';y = y(:);
        nop = isnan(x) | isnan(y);
        [rho, yfit] = LinearRegression(x(~nop),y(~nop));
        [c,p]= corr(x(~nop),y(~nop));
        plot(x(~nop),yfit,'c');
        ylim([0 3]);
        
        
    end
end
%% Plot difference performance light-off light-on, 10-90, ctl grouped VC and AC

maxday = 21;
days = 1:maxday;
savefig = false;
AC = 1;VC =2;
for i=3
    ok = explist==i;
    fig = figure;
    set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
    orient(fig,'landscape');
    taken = (status==1)';        
    groups = genlist(ok & taken);    
    lists = struclist(ok & taken);    

    subplot(2,2,1);hold on;
    data = opc(ok & taken,days)-rpc(ok & taken,days);  
%     data = opc(ok & taken,days);  
    plot(nanmean(data(groups==2 & lists==AC,:))','k','linewidth',3);
    plot(nanmean(data(groups==2 & lists==AC,:))'+nansem(data(groups==2,:))','k');
    plot(nanmean(data(groups==2 & lists==AC,:))'-nansem(data(groups==2,:))','k');
    plot(nanmean(data(groups==1 & lists==AC,:))','c','linewidth',3);
    plot(nanmean(data(groups==1 & lists==AC,:))'+nansem(data(groups==1 & lists==AC,:))','c');
    plot(nanmean(data(groups==1 & lists==AC,:))'-nansem(data(groups==1 & lists==AC,:))','c');
    PlotHVLines(0,'h','k');
    ylim([-50 10]);
%     ylim([50 100]);
    ylabel('Diff percent correct');
    xlabel([0 22]);xlabel('Days');
    title([expNames{i} '-' strucNames{1} ' (ctl=' num2str(sum(groups==2 & lists==AC)) ...
        ', test=' num2str(sum(groups==1 & lists==AC)) ') - AC only']);   
    d = [data(groups==2 & lists==AC,:)' data(groups==1 & lists==AC,:)'];
    g = [ones(sum(groups==2 & lists==AC),1);2*ones(sum(groups==1 & lists==AC),1)]';
    g = repmat(g,size(d,1),1);
%     days = repmat((1:21)',1,size(d,2));        
%     [p,~,stats] = anovan(d(:),{g(:) days(:)});
%         multcompare(stats);
%     [p,~,stats] = anova1(d(:),g(:),'off');
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
%         [~,tablep(pp)] = ttest2(data(groups==2 & lists==AC,pp)',data(groups==1 & lists==AC,pp)');
        tablep(pp) = rndttest(data(groups==2 & lists==AC,pp)',data(groups==1 & lists==AC,pp)');
    end
    if sum(tablep<0.05)~=0
        plot(find(tablep<0.05),-20,'r*');
    end
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
        [~,tablep(pp)] = ttest(data(groups==2 & lists==AC,pp)');
    end
    if sum(tablep<0.05)~=0
        plot(find(tablep<0.05),-25,'k*');
    end
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
        [~,tablep(pp)] = ttest(data(groups==1 & lists==AC,pp)');
    end
    if sum(tablep<0.05)~=0
        plot(find(tablep<0.05),-30,'c*');
    end
    
    subplot(2,2,2);hold on;
    plot(nanmean(data(groups==2,:))','k','linewidth',3);
    plot(nanmean(data(groups==2,:))'+nansem(data(groups==2,:))','k');
    plot(nanmean(data(groups==2,:))'-nansem(data(groups==2,:))','k');
    plot(nanmean(data(groups==1 & lists==AC,:))','c','linewidth',3);
    plot(nanmean(data(groups==1 & lists==AC,:))'+nansem(data(groups==1 & lists==AC,:))','c');
    plot(nanmean(data(groups==1 & lists==AC,:))'-nansem(data(groups==1 & lists==AC,:))','c');
    PlotHVLines(0,'h','k');
    ylim([-50 10]);
%     ylim([50 100]);
    ylabel('Diff percent correct');
    xlabel([0 22]);xlabel('Days');
    title([expNames{i} '-' strucNames{1} ' (ctl=' num2str(sum(groups==2)) ...
        ', test=' num2str(sum(groups==1 & lists==AC)) ') - AC + VC ctl ']);   
    d = [data(groups==2,:)' data(groups==1 & lists==AC,:)'];
    g = [ones(sum(groups==2),1);2*ones(sum(groups==1 & lists==AC),1)]';
    g = repmat(g,size(d,1),1);
%     days = repmat((1:21)',1,size(d,2));        
%     [p,~,stats] = anovan(d(:),{g(:) days(:)});
%         multcompare(stats);
%     [p,~,stats] = anova1(d(:),g(:),'off');
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
%         [~,tablep(pp)] = ttest2(data(groups==2,pp)',data(groups==1 & lists==AC,pp)');
        tablep(pp) = rndttest(data(groups==2,pp)',data(groups==1 & lists==AC,pp)');
    end
    plot(find(tablep<0.05),-20,'r*');
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
        [~,tablep(pp)] = ttest(data(groups==2,pp)');
    end
    plot(find(tablep<0.05),-25,'k*');
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
        [~,tablep(pp)] = ttest(data(groups==1 & lists==AC,pp)');
    end
    plot(find(tablep<0.05),-30,'c*');
    
    subplot(2,2,3);hold on;
    plot(nanmean(data(groups==2 | (groups==1 & lists==VC),:))','k','linewidth',3);
    plot(nanmean(data(groups==2 | (groups==1 & lists==VC),:))'+nansem(data(groups==2,:))','k');
    plot(nanmean(data(groups==2 | (groups==1 & lists==VC),:))'-nansem(data(groups==2,:))','k');
    plot(nanmean(data(groups==1 & lists==AC,:))','c','linewidth',3);
    plot(nanmean(data(groups==1 & lists==AC,:))'+nansem(data(groups==1 & lists==AC,:))','c');
    plot(nanmean(data(groups==1 & lists==AC,:))'-nansem(data(groups==1 & lists==AC,:))','c');
    PlotHVLines(0,'h','k');
    ylim([-50 10]);
%     ylim([50 100]);
    ylabel('Diff percent correct');
    xlabel([0 22]);xlabel('Days');
    title([expNames{i} '-' strucNames{1} ' (ctl=' num2str(sum(groups==2 | (groups==1 & lists==VC))) ...
        ', test=' num2str(sum(groups==1 & lists==AC)) ') - AC + VC ctl&test']);   
    d = [data(groups==2 | (groups==1 & lists==VC),:)' data(groups==1 & lists==AC,:)'];
    g = [ones(sum(groups==2 | (groups==1 & lists==VC)),1);2*ones(sum(groups==1 & lists==AC),1)]';
    g = repmat(g,size(d,1),1);
%     days = repmat((1:21)',1,size(d,2));        
%     [p,~,stats] = anovan(d(:),{g(:) days(:)});
%         multcompare(stats);
%     [p,~,stats] = anova1(d(:),g(:),'off');
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
%         [~,tablep(pp)] = ttest2(data(groups==2 | (groups==1 & lists==VC),pp)',...
%             data(groups==1 & lists==AC,pp)');
        tablep(pp) = rndttest(data(groups==2 | (groups==1 & lists==VC),pp)',...
            data(groups==1 & lists==AC,pp)');
    end
    plot(find(tablep<0.05),-20,'r*');
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
        [~,tablep(pp)] = ttest(data(groups==2 | (groups==1 & lists==VC),pp)');
    end
    plot(find(tablep<0.05),-25,'k*');
    tablep = nan(size(d,1),1);
    for pp=1:size(d,1)
        [~,tablep(pp)] = ttest(data(groups==1 & lists==AC,pp)');
    end
    plot(find(tablep<0.05),-30,'c*');
    
    if savefig
        cd(pathsave);
        t = [expNames{i} '-' strucNames{1} '-10-90-VCandAC_ctl_grouped'];
        saveas(fig,[t '.pdf']);
        close(fig);
    end  
end

%% Plot cumulative difference FA rate light-off light-on, 10-90, ctl grouped VC and AC

maxday = 16;
days = 1:maxday;
savefig = false;
AC = 1;VC =2;
for i=2
    ok = explist==i;
    fig = figure;
    set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
    orient(fig,'landscape');
    taken = (status==1)';        
    groups = genlist(ok & taken);    
    lists = struclist(ok & taken);    
    
    data = ofa(ok & taken,days)-rfa(ok & taken,days);
    datahit = rhit(ok & taken,days);
    
    [~,findmaxprobe] = max(ppc,[],2);
    
    % POST ACQUI
    subplot(4,3,1);hold on;
    dctl = data(groups==2,:)';
    start_ctl = findmaxprobe(groups==2);
    cumsumdctl = nan(size(dctl));
    for cc=1:size(dctl,2)
        le = size(dctl,1)-start_ctl(cc)+1;
        cumsumdctl(1:le,cc) = cumsum(dctl(start_ctl(cc):end,cc));
    end
    plot(cumsumdctl);
    ylim([-3 3]);    
    xlabel('Days post acquisition');
    title([expNames{i} ' - AC+VC ctl (n=' num2str(sum(groups==2)) ')']);
    ylabel('Cum sum of FA rate difference ON-OFF');
    
    subplot(4,3,2);hold on;
    dtestAC = data(groups==1 & lists==AC,:)';
    start_test = findmaxprobe(groups==1 & lists==AC);
    cumsumdtest = nan(size(dtestAC));
    for cc=1:size(dtestAC,2)
        le = size(dtestAC,1)-start_test(cc)+1;
        cumsumdtest(1:le,cc) = cumsum(dtestAC(start_test(cc):end,cc));
    end
    plot(cumsumdtest);
    ylim([-3 3]);    
    xlabel('Days post acquisition');
    title([expNames{i} ' - AC test (n=' num2str(sum(groups==1 & lists==AC)) ')']);
    
    subplot(4,3,3);hold on;
    plot(nanmean(cumsumdctl,2),'k','linewidth',2);
    plot(nanmean(cumsumdctl,2)+nansem(cumsumdctl')','k');
    plot(nanmean(cumsumdctl,2)-nansem(cumsumdctl')','k');
    plot(nanmean(cumsumdtest,2),'c','linewidth',2);
    plot(nanmean(cumsumdtest,2)+nansem(cumsumdtest')','c');
    plot(nanmean(cumsumdtest,2)-nansem(cumsumdtest')','c');
    ylim([-3 3]);    
    xlabel('Days post acquisition');
    
    % DAYS
    subplot(4,3,4);hold on;
    dctl = data(groups==2,:)';
    cumsumdctl = cumsum(dctl);
    plot(cumsumdctl);
    ylim([-3 3]);    
    xlabel('Days');
    title([expNames{i} ' - AC+VC ctl (n=' num2str(sum(groups==2)) ')']);
    ylabel('Cum sum of FA rate difference ON-OFF');
    
    subplot(4,3,5);hold on;
    dtestAC = data(groups==1 & lists==AC,:)';
    cumsumdtest = cumsum(dtestAC);
    plot(cumsumdtest);
    ylim([-3 3]);    
    xlabel('Days');
    title([expNames{i} ' - AC test (n=' num2str(sum(groups==1 & lists==AC)) ')']);
    
    subplot(4,3,6);hold on;
    plot(nanmean(cumsumdctl,2),'k','linewidth',2);
    plot(nanmean(cumsumdctl,2)+nansem(cumsumdctl')','k');
    plot(nanmean(cumsumdctl,2)-nansem(cumsumdctl')','k');
    plot(nanmean(cumsumdtest,2),'c','linewidth',2);
    plot(nanmean(cumsumdtest,2)+nansem(cumsumdtest')','c');
    plot(nanmean(cumsumdtest,2)-nansem(cumsumdtest')','c');
    ylim([-3 3]);    
    xlabel('Days');
    
    % POST ACQUI, HIGH HIT ONLY (>0.95)
    subplot(4,3,7);hold on;
    dctl = data(groups==2,:)';
    start_ctl = findmaxprobe(groups==2);
    highHIT_ctl = datahit(groups==2,:)'>=0.95;    
    cumsumdctl = nan(size(dctl));
    for cc=1:size(dctl,2)
        ddd = dctl(start_ctl(cc):end,cc);
        okk = highHIT_ctl(start_ctl(cc):end,cc);
        cumsumdctl(1:sum(okk),cc) = cumsum(ddd(okk));
    end
    plot(cumsumdctl);
    ylim([-3 3]);    
    xlabel('Days post acquisition');
    title([expNames{i} ' - AC+VC ctl (n=' num2str(sum(groups==2)) ')']);
    ylabel('Cum sum of FA diff, HIT>=0.95');
    
    subplot(4,3,8);hold on;
    dtestAC = data(groups==1 & lists==AC,:)';
    start_test = findmaxprobe(groups==1 & lists==AC);
    highHIT_test = datahit(groups==1 & lists==AC,:)'>=0.95;    
    cumsumdtest = nan(size(dtestAC));
    for cc=1:size(dtestAC,2)
        ddd = dtestAC(start_test(cc):end,cc);
        okk = highHIT_test(start_test(cc):end,cc);
        cumsumdtest(1:sum(okk),cc) = cumsum(ddd(okk));
    end
    plot(cumsumdtest);
    ylim([-3 3]);    
    xlabel('Days post acquisition');
    title([expNames{i} ' - AC test (n=' num2str(sum(groups==1 & lists==AC)) ')']);
    
    subplot(4,3,9);hold on;
    plot(nanmean(cumsumdctl,2),'k','linewidth',2);
    plot(nanmean(cumsumdctl,2)+nansem(cumsumdctl')','k');
    plot(nanmean(cumsumdctl,2)-nansem(cumsumdctl')','k');
    plot(nanmean(cumsumdtest,2),'c','linewidth',2);
    plot(nanmean(cumsumdtest,2)+nansem(cumsumdtest')','c');
    plot(nanmean(cumsumdtest,2)-nansem(cumsumdtest')','c');
    ylim([-3 3]);    
    xlabel('Days post acquisition');
    
    % DAYS
    subplot(4,3,10);hold on;
    plot(data(groups==2,:)');
    ylim([-1 0.5]);    
    xlabel('Days');
    title([expNames{i} ' - AC+VC ctl (n=' num2str(sum(groups==2)) ')']);
    ylabel('FA rate difference ON-OFF');
    
    subplot(4,3,11);hold on;
    plot(data(groups==1 & lists==AC,:)');
    ylim([-1 0.5]);    
    xlabel('Days');
    title([expNames{i} ' - AC test (n=' num2str(sum(groups==1 & lists==AC)) ')']);
    
    subplot(4,3,12);hold on;
    plot(nanmean(data(groups==2,:)),'k','linewidth',2);
    plot(nanmean(data(groups==2,:))+nansem(data(groups==2,:)),'k');
    plot(nanmean(data(groups==2,:))-nansem(data(groups==2,:)),'k');
    plot(nanmean(data(groups==1 & lists==AC,:)),'c','linewidth',2);
    plot(nanmean(data(groups==1 & lists==AC,:))+nansem(data(groups==1 & lists==AC,:)),'c');
    plot(nanmean(data(groups==1 & lists==AC,:))-nansem(data(groups==1 & lists==AC,:)),'c');
    ylim([-1 0.5]);    
    xlabel('Days');
 
    if savefig
        cd(pathsave);
        t = [expNames{i} '-DiffFA'];
        saveas(fig,[t '.pdf']);
        close(fig);
    end  
end

%% Comparison between controls in all experiments

VAR = {rhit ohit phit; rfa ofa pfa; rdprime odprime pdprime; rpc opc ppc;...
    rcriterion ocriterion pcriterion};
ylims = [[0 1];[0 1];[-1 5];[40 100];[-2 2]];
nExp = numel(unique(explist));
condNames = {'Reinf OFF','Reinf ON','Probe'};
varNames = {'hit','fa','dprime','PercentCorrect','Criterion'};
savefig = true;

fig = figure;
set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
orient(fig,'landscape')
colors = {'k','g','r','b'};
styles = {'-',':'};
names = {};
for i=1:nExp
    ok = explist==i+1;
    nStruc = numel(unique(struclist(ok)));
    for j=1:nStruc
        ok = explist==i & struclist==j;
        for v=1:size(VAR,1)             
            
            var = VAR(v,:);
            for c=1:3 % loop trough conditions, i.e. reinf, opto, probe
                if c==3 
                    taken = (status==1 | status==2)';
                else
                    taken = (status==1)';
                end
                groups = genlist(ok & taken);            
                data = var{c}(ok & taken,:);
                subplot(3,5,v+((c-1)*size(VAR,1)));hold on;
                subplot(2,size(var,2),c);hold on;
                plot(nanmean(data(groups==2,:))','color',colors{i},'linestyle',styles{j},'linewidth',1);
                plot(nanmean(data(groups==2,:))'+nansem(data(groups==2,:))','color',colors{i});
                plot(nanmean(data(groups==2,:))'-nansem(data(groups==2,:))','color',colors{i});
                
                title([expNames{i} '-' strucNames{j} ', ' condNames{c} ' (ctl=' num2str(sum(groups==2))...
                    ', test=' num2str(sum(groups==1)) ')']);
                if v==1
                    ylabel([condNames{c} ' - ' varNames{v}]);
                else
                    ylabel(varNames{v});
                end
                xlabel('Days');
                ylim(ylims(v,:));
                xlim([0 30]);
%                 if v==5, PlotHVLines(0,'h','k'); end
            end
        end   
        names = [names [expNames{i} '-' strucNames{j} ' (' num2str(sum(groups==2)) ')']];
    end
end
legend(names)
if savefig
    cd(pathsave);
    t = 'Comparison-Controls-All-Exp';
    saveas(fig,[t '.pdf']);
    close(fig);
end 

% % stats
% d = VAR{4,1};
% ok = status==1;
% code = [explist' struclist'];
% [~,w] = unique(code,'rows');
% newcode = nan(size(code,1),1);
% for l=1:length(w)
%    newcode(w(l):end) = l*ones(abs(diff([size(code,1);w(l)]))+1,1);
% end
% newcode = repmat(newcode,1,size(d,2));
% dd = d(ok & genlist'==2,:);
% gg = newcode(ok & genlist'==2,:);
% [p,~,stats] = anova1(dd(:),gg(:));
% multcompare(stats)


i=2; % 90-10 experiment
ok = explist==i;
maxday = 21;
criteriamethod = 'percent'; % 'percent' or 'dprime' or 'maxpercent'
criteriapercent = 70;
criteriadprime = 1.5;
nStruc = numel(unique(struclist(ok)));
fig = figure;
% set(fig, 'units','normalized','outerposition',[0 0.5 0.5 0.5]);
% orient(fig,'landscape');
for j=1:nStruc
    ok = explist==i & struclist==j;
    taken = (status==1)';
    groups = genlist(ok & taken);
    n = sum(ok & taken);
    if strcmp(criteriamethod,'percent')
        v = 4; % percent correct
        var = VAR(v,:);
        thresh = var{3}(ok & taken,1:maxday)>=criteriapercent;
    elseif strcmp(criteriamethod,'dprime')
        v = 3; % d'
        var = VAR(v,:);
        thresh = var{3}(ok & taken,1:maxday)>=criteriadprime;
    elseif strcmp(criteriamethod,'maxpercent')
        v = 4; % percent correct
        var = VAR(v,:);
        [~,where] = max(var{3}(ok & taken,1:maxday),[],2);
        thresh = zeros(n,maxday);
        for t=1:n
            thresh(t,where(t):end)=1;
        end
        thresh = logical(thresh);
    end
    day_criteria = nan(n,1);
    for t=1:n
        dc = strfind(thresh(t,:),[0 1]);
        if ~isempty(dc)
            day_criteria(t) = dc(1);
        elseif strfind(thresh(t,:),1)==1
            day_criteria(t) = 1;
        end
    end
    dc_ctl = day_criteria(groups==2);
    dc_test = day_criteria(groups==1);
    subplot(2,2,1);hold on;
    histogram(dc_ctl,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','k');
    histogram(dc_test,'Normalization','cdf','DisplayStyle','stairs','EdgeColor','c');
    ylabel('Cum fraction of mice');
    if strcmp(criteriamethod,'percent')
        xlabel(['Probe day > ' num2str(criteriapercent) '%']);
        title([expNames{i} '-' strucNames{j} ' (ctl=' num2str(length(dc_ctl)) ...
            ', test=' num2str(length(dc_test)) ')']);
    elseif strcmp(criteriamethod,'dprime')
        xlabel(['Probe day d" > ' num2str(criteriadprime)]);
        title([expNames{i} '-' strucNames{j} ' (ctl=' num2str(length(dc_ctl)) ...
            ', test=' num2str(length(dc_test)) ')']);
    elseif strcmp(criteriamethod,'maxpercent')
        xlabel('Probe day max %');
        title([expNames{i} '-' strucNames{j} ' (ctl=' num2str(length(dc_ctl)) ...
            ', test=' num2str(length(dc_test)) ')']);
    end
    subplot(2,2,2);hold on;
    v = 4; % percent correct
    var = VAR(v,:);
    percentReinflightOFF = var{1}(ok & taken,1:maxday);
    align_perf = nan(n,maxday);
    for t=1:n
        si = maxday-day_criteria(t)+1;
        if isnan(si), continue, end
        align_perf(t,1:si) = percentReinflightOFF(t,day_criteria(t):end);
    end
    if j == 1
        plot(align_perf(groups==2,:)','k');
        plot(align_perf(groups==1,:)','c');
        ylim([40 100]);ylabel('% correct');
        title('Reinforce light-off');
        if strcmp(criteriamethod,'percent')
            xlabel(['Day from probe day > ' num2str(criteriapercent) '%']);
        elseif strcmp(criteriamethod,'dprime')
            xlabel(['Day from probe day d" > ' num2str(criteriadprime)]);
        elseif strcmp(criteriamethod,'maxpercent')
            xlabel('Day from max % probe');
        end
        subplot(2,2,3);hold on;
        plot(nanmean(align_perf(groups==2,:)),'k','linewidth',2);
        plot(nanmean(align_perf(groups==2,:))+nansem(align_perf(groups==2,:)),'k');
        plot(nanmean(align_perf(groups==2,:))-nansem(align_perf(groups==2,:)),'k');
        plot(nanmean(align_perf(groups==1,:)),'Color', [0.5 0.5 0.5],'linewidth',2);
        plot(nanmean(align_perf(groups==1,:))+nansem(align_perf(groups==1,:)),'c');
        plot(nanmean(align_perf(groups==1,:))-nansem(align_perf(groups==1,:)),'c');
        ylim([40 100]);ylabel('% correct');
        title('Reinforce light-off');
        if strcmp(criteriamethod,'percent')
            xlabel(['Day from probe day > ' num2str(criteriapercent) '%']);
        elseif strcmp(criteriamethod,'dprime')
            xlabel(['Day from probe day d" > ' num2str(criteriadprime)]);
        elseif strcmp(criteriamethod,'maxpercent')
            xlabel('Day from max % probe');
        end
    else
        plot(nanmean(align_perf(groups==2,:)),'k','linewidth',2);
    plot(nanmean(align_perf(groups==2,:))+nansem(align_perf(groups==2,:)),'k');
    plot(nanmean(align_perf(groups==2,:))-nansem(align_perf(groups==2,:)),'k');
    plot(nanmean(align_perf(groups==1,:)),'Color', [0.5 0.5 0.5],'linewidth',2);
    plot(nanmean(align_perf(groups==1,:))+nansem(align_perf(groups==1,:)),'c');
    plot(nanmean(align_perf(groups==1,:))-nansem(align_perf(groups==1,:)),'c');
    ylim([40 100]);ylabel('% correct');
    title('Reinforce light-off');
    if strcmp(criteriamethod,'percent')
        xlabel(['Day from probe day > ' num2str(criteriapercent) '%']);
    elseif strcmp(criteriamethod,'dprime')
        xlabel(['Day from probe day d" > ' num2str(criteriadprime)]);
	elseif strcmp(criteriamethod,'maxpercent')
        xlabel('Day from max % probe');
    end
    end
end


maxdays = 30;
rhit = nan(nSubj,maxdays);
rfa = nan(nSubj,maxdays);
phit = nan(nSubj,maxdays);
pfa = nan(nSubj,maxdays);
ohit = nan(nSubj,maxdays);
ofa = nan(nSubj,maxdays);
cofffa = nan(nSubj,maxdays);
confa = nan(nSubj,maxdays);
for s=1:nSubj
    if isempty(tempsubjrates{s}), continue, end
    n = size(tempsubjrates{s},1);
    rhit(s,1:n) = tempsubjrates{s}(:,1);
    rfa(s,1:n) = tempsubjrates{s}(:,2);
    phit(s,1:n) = tempsubjrates{s}(:,3);
    pfa(s,1:n) = tempsubjrates{s}(:,4);
    ohit(s,1:n) = tempsubjrates{s}(:,5);
    ofa(s,1:n) = tempsubjrates{s}(:,6);
    cofffa(s,1:n) = tempsubjrates{s}(:,9);
    confa(s,1:n) = tempsubjrates{s}(:,10);
end
% Compute percent correct
rpc = (rhit+(1-rfa))/2*100;
opc = (ohit+(1-ofa))/2*100;
starts = find(sum(isnan(opc))~=nSubj);
starts(end)=[]; % remove light outside the brain
figure;
for s=1:nSubj
    subplot(2,2,s);hold on;
    plot([rpc(s,starts)' opc(s,starts)']','k.-');
    xlim([0.5 2.5]);
    set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
    ylim([0 100]);ylabel('Proportion correct');
    [~,p] = ttest(rpc(s,starts)',opc(s,starts)');
    title([subjlist{s} ', p=' num2str(p)]);
end
figure;hold on;
subplot(2,2,1);hold on;
lightoff = rpc(:,starts)';
lighton = opc(:,starts)';
plot([lightoff(:) lighton(:)]','k.-');
plot([1 2],[nanmean(lightoff(:)) nanmean(lighton(:))],'b.','markersize',15);
xlim([0.5 2.5]);
set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
ylim([0 100]);ylabel('Proportion correct');
subplot(2,2,2);hold on;
colors = {'c','y','b','g'};
plot([1 2],[nanmean(lightoff(:)) nanmean(lighton(:))],'k.','markersize',15);
for s=1:nSubj
    lightoff = nanmean(rpc(s,starts),2);
    lighton = nanmean(opc(s,starts),2);
    plot([lightoff lighton]','.-','color',colors{s});
end
xlim([0.5 2.5]);
set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
ylim([0 100]);ylabel('Proportion correct');
subplot(2,2,3);hold on;
title('Light outside the brain');
outside = find(sum(isnan(opc))~=nSubj);
outside = outside(end);
lightoff = rpc(:,outside)';
lighton = opc(:,outside)';
plot([lightoff(:) lighton(:)]','k.-');
plot([1 2],[nanmean(lightoff(:)) nanmean(lighton(:))],'b.','markersize',15);
xlim([0.5 2.5]);
set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
ylim([0 100]);ylabel('Proportion correct');