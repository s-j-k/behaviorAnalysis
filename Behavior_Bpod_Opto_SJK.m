function Behavior_Bpod_Opto_SJK(cohort)

%MGB IC opto cohort

switch cohort
    case 1
    subjlist={'sk158','sk160','sk161', 'sk162','sk163','sk164'}; %PV Cre
    explist=[2 2 2 2 1 1]; % 2 is ctl
    pathsave='O:\sjk\Behavior\optoMGBIC_1\';
    
    case 2
    subjlist={'sk170','sk175'}; %GTACR
    explist=[2 1];
    pathsave='O:\sjk\Behavior\OptoMGBIC_2\';

    case 3
    pathsave='O:\sjk\Behavior\optoMGBIC_3\';
    subjlist={'sk176','sk177','sk178','sk179'}; %GTACR
    explist=[1 2 1 2];
    
    case 4
    pathsave='O:\sjk\Behavior\MGBIC_4\';
    subjlist={'sk182','sk183','sk184','sk185','sk186','sk187','sk188','sk189'}; %GTACR
    explist=[1 2 1 1 1 1 1 2];
    
    case 5
    pathsave='O:\sjk\Behavior\MGBIC_5\';
%     pathsave='C:\Users\sjkim1\Desktop\OptoData\MGBIC_5\';
    subjlist={'sk190','sk191','sk192','sk193','sk194','sk195','sk196','sk197'}; %GTACR
%     explist=[2 2 1 1 1 1 1 1];
    explist=[1 1 1 1 1 2 1 1];
    
    case 6
    pathsave='O:\sjk\Behavior\MGBIC_6\';
%     pathsave='C:\Users\sjkim1\Desktop\OptoData\MGBIC_5\';
    subjlist={'sk198','sk199','sk200','sk201','sk202','sk203','sk204','sk205'}; %GTACR
    explist=[1 1 1 1 2 1 1 1];    
   case 7
    pathsave='O:\sjk\Behavior\IC_cohort_1\';
    subjlist={'sk218','sk219','sk220'}; %GTACR
    explist=[1 1 1];    

    case 8
    pathsave='O:\sjk\Behavior\IC_cohort_2\';
%     pathsave='C:\Users\sjkim1\Desktop\OptoData\MGBIC_5\';
    subjlist={'sk223','sk224','sk225','sk226'}; %GTACR
    explist=[1 1 1 1]
    
    otherwise
        disp('Cohort not found');    
end

% subjlist = {'SK49','SK50','SK51','SK52','SK53','SK54','SK55','SK56'};
% explist = [1 1 1 1 1 2 2 2]; % 1 = test, 2 = ctl
% pathsave='Z:\su\DATA\behaviorData\opto\cohort2\';

expnames = {'TEST', 'CTL'};
path=pathsave;

nSubj = length(subjlist);
plot_indiv_data = true;
plot_lick_psth = true;
savefig = true;
filetype='fig';

nFiles = 17; % = nb of days of behavior

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
optomeanMat={};ii=2;
optomeanMat{1,ii}='RHit'; optomeanMat{1,ii+1}='RFA';
optomeanMat{1,ii+2}='Full Opto Hit'; optomeanMat{1,ii+3}='Full Opto FA';
optomeanMat{1,ii+4}='Tone Opto Hit'; optomeanMat{1,ii+5}='Tone Opto FA';
optomeanMat{1,ii+6}='Choice Opto Hit'; optomeanMat{1,ii+7}='Choice Opto FA';
optomeanMat{1,ii+8}= 'MGB R Hit';optomeanMat{1,ii+9}='MGB R FA';
optomeanMat{1,ii+10}='MGB Full O Hit'; optomeanMat{1,ii+11}='MGB Full O FA';
optomeanMat{1,ii+12}='MGB Tone Opto Hit';  optomeanMat{1,ii+13}='MGB Tone Opto FA';
optomeanMat{1,ii+14}='MGB Choice Opto Hit'; optomeanMat{1,ii+15}='MGB Choice Opto FA';
optomeanMat{1,ii+16}= 'IC R Hit'; optomeanMat{1,ii+17}='IC R FA';
optomeanMat{1,ii+18}='IC Full Opto Hit';optomeanMat{1,ii+19}='IC Full Opto FA';
optomeanMat{1,ii+20}='IC Tone Opto Hit';  optomeanMat{1,ii+21}='IC Tone Opto FA';
optomeanMat{1,ii+22}='IC Choice Opto Hit'; optomeanMat{1,ii+23}='IC Choice Opto FA';
optomeanMat{1,ii+24}='Matrix Variable';optomeanMat{1,ii+25}='Rates Variable';

for nbsubj = 1:nSubj % through animals
    % Localize data
    subj = subjlist{nbsubj}; 
  
    if str2double(subj(end-2:end))>100
        subj = lower(subj);
    end  
%     if explist(nbsubj)>= 4 || nbsubj > 65
%         subjPath = [path subj '\GNG_LaserSineWavePhasicTone_Celine_ZZAW\Session Data\']; 
%     elseif str2double(subj(end-2:end)) >= 123
%         subjPath = [path subj '\GNG_LaserSineWavePhasicTone_Celine_ZZAW\Session Data\']; 
%     else
%         subjPath = [path subj '\GNG_LaserSineWavePhasicTone_Celine\Session Data\'];
%     end
%     subjPath = [path subj '\GNG_LaserSineWavePhasicTone_Celine_ZZAW\Session Data\'];
%     subjPath = [path subj '\GNG_MultiProbOpto_SJK\Session Data\'];
    subjPath = [path subj '\GNG_MultiProbOptoChoice_SJK\Session Data\'];
%     subjPath = [path subj '\LickingGNG\Session Data\'];

    files = dir([subjPath '*.mat']); 
    nFiles = length(files);
    if nFiles == 0; end
    
    filenames = cell(nFiles,1); 
    licks = cell(nFiles,1);
    matrix = []; 
    SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
    START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
    toremove = 0;
    for i=1:nFiles

        if ~exist([subjPath 'bad']) % TRYING TO remove files that are too
%         small so I don't do it manually...
            mkdir(subjPath, 'bad') %makes directory to add here files that have less than 10 trials to be excluded from analysis
        end 
        cd(subjPath)
        matfiles2 = dir('*.mat'); %finds all GNG sessions for the selected animal

        %to make sure files with less than 270 trials are not taken into
        %consideration
        for maf=1:length(matfiles2)
            clear namemat
            namemat=matfiles2(maf).name;
            load(namemat);
            if SessionData.nTrials<270
                movefolder= ([subjPath 'bad']);
                movefile(namemat, movefolder) %moves files to sessionOut folder per animal
            end
            cd(subjPath);
        end
    end
    
    files = dir([subjPath '*.mat']); % do it again to get the new file list
    nFiles = length(files);
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
                if isfield(SessionData.RawEvents.Trial{1,j}.Events,'Port1In') % for bottom box, where licks were not detected
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
            elseif tempMat(j,CTXT) == 5 % tone
                if ~isnan(SessionData.RawEvents.Trial{1,j}.States.OpenValve)
                    tempMat(j,OUTCOME) = 1; tempMat(j,TONE) = 1;
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.Miss)
                    tempMat(j,OUTCOME) = 2; tempMat(j,TONE) = 1; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.Punish)
                    tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 2; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.CorrectReject) 
                    tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 2;
                end
            elseif tempMat(j,CTXT) == 6 % choice
                if ~isnan(SessionData.RawEvents.Trial{1,j}.States.OpenValve)
                    tempMat(j,OUTCOME) = 1; tempMat(j,TONE) = 1;
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.Miss)
                    tempMat(j,OUTCOME) = 2; tempMat(j,TONE) = 1; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.PunishFirst)
                    tempMat(j,OUTCOME) = 3; tempMat(j,TONE) = 2; 
                elseif ~isnan(SessionData.RawEvents.Trial{1,j}.States.CorrectReject) 
                    tempMat(j,OUTCOME) = 4; tempMat(j,TONE) = 2;
                end
                
            end
        end 
        licks{numSess} = [licks{numSess};templicks];
        matrix = [matrix;tempMat];
    end
    conditionNum = 18; % this is 14 when there is only one type of opto condition, but 18 if there are 3 oopto conditions
    ndays = max(matrix(:,SESS)); % if this line throws an error check the protocol path, lines 120-130-ish
    rates = nan(ndays,conditionNum);
    
    nctxt = nan(ndays,conditionNum);
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
        rfhit = sum(tempMat(:,OUTCOME)==1 & (tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2 | tempMat(:,CTXT)==5 | tempMat(:,CTXT)==6) / sum((tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2 | tempMat(:,CTXT)==5 | tempMat(:,CTXT)==6) & tempMat(:,TONE)==1)); % hit in reinforced light off + on
        phit2_1 = sum(tempMat(:,OUTCOME)==1 & probeIdx1) / sum(probeIdx1 & tempMat(:,TONE)==1); % hit in probe, dividing probe bloc in 2
        phit2_2 = sum(tempMat(:,OUTCOME)==1 & probeIdx2) / sum(probeIdx2 & tempMat(:,TONE)==1); % hit in probe, dividing probe bloc in 2
        phit2 = [phit2_1 phit2_2];
        othit = sum(tempMat(:,OUTCOME)==1 & tempMat(:,CTXT)==5) / sum(tempMat(:,CTXT)==5 & tempMat(:,TONE)==1); % hit in opto tone (i.e. reinforced light on)
        ochit = sum(tempMat(:,OUTCOME)==1 & tempMat(:,CTXT)==6) / sum(tempMat(:,CTXT)==6 & tempMat(:,TONE)==1); % hit in opto choice (i.e. reinforced light on)
        
        rfa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==2) / sum(tempMat(:,CTXT)==2 & tempMat(:,TONE)==2);
        pfa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==0) / sum(tempMat(:,CTXT)==0 & tempMat(:,TONE)==2);
        ofa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==1) / sum(tempMat(:,CTXT)==1 & tempMat(:,TONE)==2);
        rffa = sum(tempMat(:,OUTCOME)==3 & (tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2)) / sum((tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2) & tempMat(:,TONE)==2);
        pfa2_1 = sum(tempMat(:,OUTCOME)==3 & probeIdx1) / sum(probeIdx1 & tempMat(:,TONE)==2);
        pfa2_2 = sum(tempMat(:,OUTCOME)==3 & probeIdx2) / sum(probeIdx2 & tempMat(:,TONE)==2);
        pfa2 = [pfa2_1 pfa2_2];
        otfa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==5) / sum(tempMat(:,CTXT)==5 & tempMat(:,TONE)==2); % hit in opto tone (i.e. reinforced light on)
        ocfa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==6) / sum(tempMat(:,CTXT)==6 & tempMat(:,TONE)==2); % hit in opto tone (i.e. reinforced light on)
        

        cofffa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==3) / sum(tempMat(:,CTXT)==3 & tempMat(:,TONE)==0); % catch light off FA
        conffa = sum(tempMat(:,OUTCOME)==3 & tempMat(:,CTXT)==4) / sum(tempMat(:,CTXT)==4 & tempMat(:,TONE)==0); % catch light on FA
        rates(i,:) = [rhit rfa phit pfa ohit ofa ...
                    rfhit rffa cofffa conffa phit2 pfa2 ...
                    othit otfa ochit ocfa];
                

        nctxt(i,:) = [sum(tempMat(:,CTXT)==2 & tempMat(:,TONE)==1) sum(tempMat(:,CTXT)==2 & tempMat(:,TONE)==2)... % numbers: e.g nb of HIT and nb of Target trials
            sum(tempMat(:,CTXT)==0 & tempMat(:,TONE)==1) sum(tempMat(:,CTXT)==0 & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)==1 & tempMat(:,TONE)==1) sum(tempMat(:,CTXT)==1 & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)~=0 & tempMat(:,TONE)==1) sum((tempMat(:,CTXT)==1 | tempMat(:,CTXT)==2) & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)==3 & tempMat(:,TONE)==0) sum(tempMat(:,CTXT)==4 & tempMat(:,TONE)==0)...
            sum(probeIdx1 & tempMat(:,TONE)==1) sum(probeIdx2 & tempMat(:,TONE)==1)...
            sum(probeIdx2 & tempMat(:,TONE)==2) sum(probeIdx2 & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)==5 & tempMat(:,TONE)==1) sum(tempMat(:,CTXT)==5 & tempMat(:,TONE)==2)...
            sum(tempMat(:,CTXT)==6 & tempMat(:,TONE)==1) sum(tempMat(:,CTXT)==6 & tempMat(:,TONE)==2)];

        reinf = tempMat(:,CTXT)==2; opto = tempMat(:,CTXT)==1; probe = tempMat(:,CTXT)==0;
        coff = tempMat(:,CTXT)==3; con = tempMat(:,CTXT)==4; topto = tempMat(:,CTXT)==5; copto = tempMat(:,CTXT)==6;
        
        rt_r = RelativeTimes(templicks,[tempMat(reinf,START)+tempMat(reinf,TONE_T)-1 tempMat(reinf,START)+tempMat(reinf,TONE_T) ...
            tempMat(reinf,START)+tempMat(reinf,TONE_T)+4],[-1 0 4]);
        rt_o = RelativeTimes(templicks,[tempMat(opto,START)+tempMat(opto,TONE_T)-1 tempMat(opto,START)+tempMat(opto,TONE_T) ...
            tempMat(opto,START)+tempMat(opto,TONE_T)+4],[-1 0 4]);
        rt_ot = RelativeTimes(templicks,[tempMat(topto,START)+tempMat(topto,TONE_T)-1 tempMat(topto,START)+tempMat(topto,TONE_T) ...
            tempMat(topto,START)+tempMat(topto,TONE_T)+4],[-1 0 4]);
        rt_oc = RelativeTimes(templicks,[tempMat(copto,START)+tempMat(copto,TONE_T)-1 tempMat(copto,START)+tempMat(copto,TONE_T) ...
            tempMat(copto,START)+tempMat(copto,TONE_T)+4],[-1 0 4]);
        lickhistr(i,:) = hist(rt_r,bins)/sum(reinf);
        lickhisto(i,:) = hist(rt_o,bins)/sum(opto);
        lickhistot(i,:) = hist(rt_ot,bins)/sum(topto);
        lickhistoc(i,:) = hist(rt_oc,bins)/sum(copto);
        hit = tempMat(:,OUTCOME)==1;fa = tempMat(:,OUTCOME)==3;
        rt_r_hit = RelativeTimes(templicks,[tempMat(reinf&hit,START)+tempMat(reinf&hit,TONE_T)-1 tempMat(reinf&hit,START)+tempMat(reinf&hit,TONE_T) tempMat(reinf&hit,START)+tempMat(reinf&hit,TONE_T)+4],[-1 0 4]);
        rt_o_hit = RelativeTimes(templicks,[tempMat(opto&hit,START)+tempMat(opto&hit,TONE_T)-1 tempMat(opto&hit,START)+tempMat(opto&hit,TONE_T) tempMat(opto&hit,START)+tempMat(opto&hit,TONE_T)+4],[-1 0 4]);
        rt_ot_hit = RelativeTimes(templicks,[tempMat(topto&hit,START)+tempMat(topto&hit,TONE_T)-1 tempMat(topto&hit,START)+tempMat(topto&hit,TONE_T) tempMat(topto&hit,START)+tempMat(topto&hit,TONE_T)+4],[-1 0 4]);
        rt_oc_hit = RelativeTimes(templicks,[tempMat(copto&hit,START)+tempMat(copto&hit,TONE_T)-1 tempMat(copto&hit,START)+tempMat(copto&hit,TONE_T) tempMat(copto&hit,START)+tempMat(copto&hit,TONE_T)+4],[-1 0 4]);
        rt_r_fa = RelativeTimes(templicks,[tempMat(reinf&fa,START)+tempMat(reinf&fa,TONE_T)-1 tempMat(reinf&fa,START)+tempMat(reinf&fa,TONE_T) tempMat(reinf&fa,START)+tempMat(reinf&fa,TONE_T)+4],[-1 0 4]);
        rt_o_fa = RelativeTimes(templicks,[tempMat(opto&fa,START)+tempMat(opto&fa,TONE_T)-1 tempMat(opto&fa,START)+tempMat(opto&fa,TONE_T) tempMat(opto&fa,START)+tempMat(opto&fa,TONE_T)+4],[-1 0 4]);
        rt_ot_fa = RelativeTimes(templicks,[tempMat(topto&fa,START)+tempMat(topto&fa,TONE_T)-1 tempMat(topto&fa,START)+tempMat(topto&fa,TONE_T) tempMat(topto&fa,START)+tempMat(topto&fa,TONE_T)+4],[-1 0 4]);
        rt_oc_fa = RelativeTimes(templicks,[tempMat(copto&fa,START)+tempMat(copto&fa,TONE_T)-1 tempMat(copto&fa,START)+tempMat(copto&fa,TONE_T) tempMat(copto&fa,START)+tempMat(copto&fa,TONE_T)+4],[-1 0 4]);
        
        rt_r_miss = RelativeTimes(templicks,[tempMat(reinf&~hit,START)+tempMat(reinf&~hit,TONE_T)-1 tempMat(reinf&~hit,START)+tempMat(reinf&~hit,TONE_T) tempMat(reinf&~hit,START)+tempMat(reinf&~hit,TONE_T)+5],[-1 0 5]);
        rt_r_cr = RelativeTimes(templicks,[tempMat(reinf&~fa,START)+tempMat(reinf&~fa,TONE_T)-1 tempMat(reinf&~fa,START)+tempMat(reinf&~fa,TONE_T) tempMat(reinf&~fa,START)+tempMat(reinf&~fa,TONE_T)+5],[-1 0 5]);
        rt_o_miss = RelativeTimes(templicks,[tempMat(opto&~hit,START)+tempMat(opto&~hit,TONE_T)-1 tempMat(opto&~hit,START)+tempMat(opto&~hit,TONE_T) tempMat(opto&~hit,START)+tempMat(opto&~hit,TONE_T)+5],[-1 0 5]);
        rt_ot_miss = RelativeTimes(templicks,[tempMat(topto&~hit,START)+tempMat(topto&~hit,TONE_T)-1 tempMat(topto&~hit,START)+tempMat(topto&~hit,TONE_T) tempMat(topto&~hit,START)+tempMat(topto&~hit,TONE_T)+5],[-1 0 5]);
        rt_oc_miss = RelativeTimes(templicks,[tempMat(copto&~hit,START)+tempMat(copto&~hit,TONE_T)-1 tempMat(copto&~hit,START)+tempMat(copto&~hit,TONE_T) tempMat(copto&~hit,START)+tempMat(copto&~hit,TONE_T)+5],[-1 0 5]);
        rt_o_cr = RelativeTimes(templicks,[tempMat(opto&~fa,START)+tempMat(opto&~fa,TONE_T)-1 tempMat(opto&~fa,START)+tempMat(opto&~fa,TONE_T) tempMat(opto&~fa,START)+tempMat(opto&~fa,TONE_T)+5],[-1 0 5]);
        rt_ot_cr = RelativeTimes(templicks,[tempMat(topto&~fa,START)+tempMat(topto&~fa,TONE_T)-1 tempMat(topto&~fa,START)+tempMat(topto&~fa,TONE_T) tempMat(topto&~fa,START)+tempMat(topto&~fa,TONE_T)+5],[-1 0 5]);
        rt_oc_cr = RelativeTimes(templicks,[tempMat(copto&~fa,START)+tempMat(copto&~fa,TONE_T)-1 tempMat(copto&~fa,START)+tempMat(copto&~fa,TONE_T) tempMat(copto&~fa,START)+tempMat(copto&~fa,TONE_T)+5],[-1 0 5]);
        
        lickhistr_hit(i,:) = hist(rt_r_hit,bins)/sum(reinf&hit);
        lickhistr_fa(i,:) = hist(rt_r_fa,bins)/sum(reinf&fa);
        lickhistr_miss(i,:) = hist(rt_r_miss,bins)/sum(reinf&~hit);
        lickhistr_cr(i,:) = hist(rt_r_cr,bins)/sum(reinf&~fa);
        
        lickhisto_hit(i,:) = hist(rt_o_hit,bins)/sum(opto&hit);
        lickhisto_fa(i,:) = hist(rt_o_fa,bins)/sum(opto&fa);
        lickhisto_miss(i,:) = hist(rt_o_miss,bins)/sum(opto&~hit);
        lickhisto_cr(i,:) = hist(rt_o_cr,bins)/sum(opto&~fa);
        
        lickhistot_hit(i,:) = hist(rt_ot_hit,bins)/sum(topto&hit);
        lickhistot_fa(i,:) = hist(rt_ot_fa,bins)/sum(topto&fa); 
        lickhistot_miss(i,:) = hist(rt_ot_miss,bins)/sum(topto&~hit);
        lickhistot_cr(i,:) = hist(rt_oc_cr,bins)/sum(topto&~fa);
        
        lickhistoc_hit(i,:) = hist(rt_oc_hit,bins)/sum(copto&hit);
        lickhistoc_fa(i,:) = hist(rt_oc_fa,bins)/sum(copto&fa);
        lickhistoc_miss(i,:) = hist(rt_oc_miss,bins)/sum(copto&~hit);
        lickhistoc_cr(i,:) = hist(rt_oc_cr,bins)/sum(copto&~fa);        
        
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
    lickhistothit{nbsubj} = lickhistot_hit;
    lickhistotfa{nbsubj} = lickhistot_fa;
    lickhistochit{nbsubj} = lickhistoc_hit;
    lickhistocfa{nbsubj} = lickhistoc_fa;
    lickhistphit{nbsubj} = lickhistp_hit;
    lickhistpfa{nbsubj} = lickhistp_fa;
    
    lickhistrmiss{nbsubj} = lickhistr_miss;
    lickhistrcr{nbsubj} = lickhistr_cr;
    lickhistomiss{nbsubj} = lickhisto_miss;
    lickhistocr{nbsubj} = lickhisto_cr;
    lickhistotmiss{nbsubj} = lickhistot_miss;
    lickhistotcr{nbsubj} = lickhistot_cr;
    lickhistochit{nbsubj} = lickhistoc_miss;
    lickhistoccr{nbsubj} = lickhistoc_cr;
    
    ratescorr = rates;
    for j=1:size(rates,2) % this is to correct for the d' in case the animal has 0 or 1 for any rate value
        ratescorr(ratescorr(:,j)==0,j) = 1./(2*nctxt(ratescorr(:,j)==0,j)); 
        ratescorr(ratescorr(:,j)==1,j) = 1-1./(2*nctxt(ratescorr(:,j)==1,j)); 
    end
    
    subjrates{nbsubj} = ratescorr;
    subjnbctxtoutcome{nbsubj} = nctxt;
    
    rdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,1)) - norminv(ratescorr(:,2)) )'; % reinforced light off
    pdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,3)) - norminv(ratescorr(:,4)) )'; % probe
    odprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,5)) - norminv(ratescorr(:,6)) )'; % opto (i.e. reinforced light on)
    otdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,15)) - norminv(ratescorr(:,16)) )'; % tone opto (i.e. reinforced light on)
    ocdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,17)) - norminv(ratescorr(:,18)) )'; % choice opto (i.e. reinforced light on)

    rfdprime(nbsubj,1:ndays) = ( norminv(ratescorr(:,7)) - norminv(ratescorr(:,8)) )'; % all reinforced (light on + off)
    
    rcriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,1)) + norminv(ratescorr(:,2)) )'/2; % reinforced light off
    ocriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,5)) + norminv(ratescorr(:,6)) )'/2; % opto (i.e. reinforced light on)
    pcriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,3)) + norminv(ratescorr(:,4)) )'/2; % probe
    otcriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,15)) + norminv(ratescorr(:,16)) )'/2; % opto (i.e. reinforced light on)
    occriterion(nbsubj,1:ndays) = -( norminv(ratescorr(:,17)) + norminv(ratescorr(:,18)) )'/2; % opto (i.e. reinforced light on)
    
    pdprime1(nbsubj,1:ndays) = ( norminv(ratescorr(:,11)) - norminv(ratescorr(:,12)) )'; % probe 10 1st trials
    pdprime2(nbsubj,1:ndays) = ( norminv(ratescorr(:,13)) - norminv(ratescorr(:,14)) )'; % probe 10 last trials
    
    for i=1:nFiles
        mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit latency
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
    expertdays=[];
    % lick latency histogram at only expert performance
    getLickLatHist(matrix,nbsubj,subjlist,expertdays)
    MAT{nbsubj,1} = matrix;
%     %%
    if plot_indiv_data
%         % Plot subj data: H and FA rates, dprimes, lick latency, lick rate
%         fig = figure(nbsubj); hold on;    
% 
%        % rates
%        % rhit rfa phit pfa ohit ofa rfhit rffa cofffa conffa phit2 pfa2];
% 
%         subplot(2,2,1); hold on; % H and FA rates
%         plot(1:ndays,rates(:,1),'-o','LineWidth',2,'color',[0.5843 0.8157 0.9882]);
%         plot(1:ndays,rates(:,2),'-o','LineWidth',2,'color',([0.9 0.6 0.6]));
%         nprob=5; % number of probe days before it reaches the peak
%         phitp=rates(:,3);
% %         phitp=phitp(1:nprob);
%         pfap=rates(:,4);
%         pfap=pfap(1:nprob);
%         plot(1:nprob,phitp,'-o','LineWidth',2,'color',[0.8 0.8 0.8]);
%         plot(1:nprob,pfap,'-o','LineWidth',2,'color',[0.8 0.8 0.8]);
% 
%         plot(1:ndays,rates(:,5),'b-o','LineWidth',2);
%         plot(1:ndays,rates(:,6),'r-o','LineWidth',2);
%         plot(1:ndays,rates(:,15),'b-x','LineWidth',2); % tone
%         plot(1:ndays,rates(:,16),'r-x','LineWidth',2); % tone
%         plot(1:ndays,rates(:,17),'b-^','LineWidth',2); % choice
%         plot(1:ndays,rates(:,18),'r-^','LineWidth',2); % choice
%         title([subjlist{nbsubj} ' (' expnames{explist(nbsubj)} ') - H and FA rates']);
%         ylim([0 1]);xlim([1 ndays]);
%         ylabel('Proportion');xlabel('Days');
%         
%         legend('reinf hit','reinf fa',...
%             'probe hit','probe fa','full opto hit','full opto fa',...
%             'tone opto hit','tone opto fa','choice opto hit','choice opto fa');
% 
%         subplot(2,2,2); hold on; % d primes
%         plot(1:ndays,rdprime(nbsubj,1:ndays),'-o','LineWidth',2,'color',[0.5843 0.8157 0.9882]);
%         plot(1:4,pdprime(nbsubj,1:4),'-o','LineWidth',2,'color',[0.8 0.8 0.8]);
%         plot(1:ndays,odprime(nbsubj,1:ndays),'b-o','LineWidth',2);  
%         plot(1:ndays,otdprime(nbsubj,1:ndays),'b-x','LineWidth',2);
%         plot(1:ndays,ocdprime(nbsubj,1:ndays),'b-^','LineWidth',2);
%         ylabel("d'");xlabel('Days');
%         ylim([-1 5]);xlim([1 ndays]);
%         PlotHVLines(0,'h','k'); legend('reinforced','probe','opto','tone opto','choice opto');
%     
%         subplot(2,4,5); hold on; % lick latency
%         shadedErrorBar(1:ndays,mr_hit(nbsubj,1:ndays),sr_hit(nbsubj,1:ndays),{'-o','LineWidth',2,'color','k'},0.5);
% %         shadedErrorBar(1:ndays,mp_hit(nbsubj,1:ndays),sp_hit(nbsubj,1:ndays),{'-o','LineWidth',2,'color',[0.8 0.8 0.8]},0.5);
%         shadedErrorBar(1:ndays,mo_hit(nbsubj,1:ndays),so_hit(nbsubj,1:ndays),{'-o','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mto_hit(nbsubj,1:ndays),sto_hit(nbsubj,1:ndays),{'-x','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mco_hit(nbsubj,1:ndays),sco_hit(nbsubj,1:ndays),{'-^','LineWidth',2,'color','b'},0.5);
%         
%         ylabel('HIT Lick latency (1st lick) (s)');legend('','reinforced','','full opto','','tone opto','','choice opto');
%         xlabel('Days');
%         
%         subplot(2,4,6); hold on; % lick latency
%         shadedErrorBar(1:ndays,mr_fa(nbsubj,1:ndays),sr_fa(nbsubj,1:ndays),{'-o','LineWidth',2,'color','k'},0.5);
% %         shadedErrorBar(1:ndays,mp_fa(nbsubj,1:ndays),sp_fa(nbsubj,1:ndays),{'-o','LineWidth',2,'color',[0.8 0.8 0.8]},0.5);
%         shadedErrorBar(1:ndays,mo_fa(nbsubj,1:ndays),so_fa(nbsubj,1:ndays),{'-o','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mto_fa(nbsubj,1:ndays),sto_fa(nbsubj,1:ndays),{'-x','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mco_fa(nbsubj,1:ndays),sco_fa(nbsubj,1:ndays),{'-^','LineWidth',2,'color','b'},0.5);
%         
%         ylabel('FA Lick latency (1st lick) (s)');legend('','reinforced','','full opto','','tone opto','','choice opto');
%         xlabel('Days');
%         
%         subplot(2,4,7); hold on; % lick rate FA
% %         shadedErrorBar(1:8,mlr_hit(nbsubj,1:8),slr_hit(nbsubj,1:8),{'-o','LineWidth',2,'color','k'},0.5);
%         shadedErrorBar(1:ndays,mlp_hit(nbsubj,1:ndays),slp_hit(nbsubj,1:ndays),{'-o','LineWidth',2,'color',[0.8 0.8 0.8]},0.5);
%         shadedErrorBar(1:ndays,mlo_hit(nbsubj,1:ndays),slo_hit(nbsubj,1:ndays),{'-o','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mlto_hit(nbsubj,1:ndays),slto_hit(nbsubj,1:ndays),{'-x','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mlco_hit(nbsubj,1:ndays),slco_hit(nbsubj,1:ndays),{'-^','LineWidth',2,'color','b'},0.5);
%         
%         ylabel('HIT Lick rate (Hz)');legend('','reinforced','','full opto','','tone opto','','choice opto');
%         xlabel('Days');
%         
%         subplot(2,4,8); hold on; % lick rate FA
% %         shadedErrorBar(1:8,mlr_fa(nbsubj,1:8),slr_fa(nbsubj,1:8),{'-o','LineWidth',2,'color','k'},0.5);
%         shadedErrorBar(1:ndays,mlp_fa(nbsubj,1:ndays),slp_fa(nbsubj,1:ndays),{'-o','LineWidth',2,'color',[0.8 0.8 0.8]},0.5);
%         shadedErrorBar(1:ndays,mlo_fa(nbsubj,1:ndays),slo_fa(nbsubj,1:ndays),{'-o','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mlto_fa(nbsubj,1:ndays),slto_fa(nbsubj,1:ndays),{'-x','LineWidth',2,'color','b'},0.5);
%         shadedErrorBar(1:ndays,mlco_fa(nbsubj,1:ndays),slco_fa(nbsubj,1:ndays),{'-^','LineWidth',2,'color','b'},0.5);
% 
%         ylabel('FA Lick rate (Hz)');legend('','reinforced','','full opto','','tone opto','','choice opto');
%         xlabel('Days'); 
%         

%         
%         if savefig
%             cd(pathsave);
%             saveas(fig,[subjlist{nbsubj} '-ActionAndPerformance.' filetype]);
%             close(fig);
%         end
        
        if plot_lick_psth
            fig=figure;hold on;
            for d=1:ndays
                subplot(8,8,d);hold on;
                plot(bins,lickhistr_hit(d,:),'k');
                plot(bins,lickhistr_fa(d,:),'k:');
                plot(bins,lickhisto_hit(d,:),'b');
                plot(bins,lickhisto_fa(d,:),'b:');
%                 PlotHVLines(0,'v','k');
%                 title([subjlist{nbsubj} ', Day ' num2str(d)]);
                xlabel('Time from tone onset (s)');
            end

            if savefig
                cd(pathsave);
                saveas(fig,[subjlist{nbsubj} '-LickPSTHreinforced-' num2str(nbins) 'bins.' filetype]);
                close(fig);
            end

            fig=figure;hold on;
            for d=1:ndays
                subplot(8,8,d);hold on;
                plot(bins,lickhistp_hit(d,:),'k');
                plot(bins,lickhistp_fa(d,:),'k:');
                PlotHVLines(0,'v','k');
                title([subjlist{nbsubj} ', Probe day ' num2str(d)]);
                xlabel('Time from tone onset (s)');
            end

            if savefig
                cd(pathsave);
                saveas(fig,[subjlist{nbsubj} '-LickPSTHprobe-' num2str(nbins) 'bins.' filetype]);
                close(fig);
            end
        end
        
        drawnow;
    end  
%     
    
%% make plots summarizing performance across training
    filetype='.fig';
    makesumperfplot=1;
    
    if makesumperfplot==1
        perfFig=figure;
        perfFig.Position(3:4)=[400,700];
        days=size(rates,1);
        subplot(3,1,1);hold on;
        hitcolor=[0.47,0.67,0.3];
        facolor=[0.93,0.69,0.3];
        plot(1:days,rates(1:days,1),'LineWidth',2,'color',hitcolor);
        plot(1:days,rates(1:days,2),'LineWidth',2,'color',facolor);
        legend('hit','fa','location','best');
        curtick = get(gca, 'xTick');
        xticks(unique(round(curtick)));
        ylim([0 1]);
        hitColorPlot=repmat(hitcolor,days,1);
        faColorPlot=repmat(facolor,days,1);
        title([subjlist{nbsubj} ' (' expnames{explist(nbsubj)} ') - H and FA rates']);
         ylabel('rate');xlabel('days');
         
        subplot(3,1,2);
        plot(1:days,rdprime(nbsubj,1:days),'LineWidth',2,'color',[0 0 0]); hold on;
        plot(1:days,pdprime(nbsubj,1:days),'LineWidth',2,'color',[0.7 0.7 0.7]);
        legend('Reinforced','Probe','location','best');ylim([0 5]);
        title('d prime');xlabel('days');
        curtick = get(gca, 'xTick');
        xticks(unique(round(curtick)));
        subplot(3,1,3);qqq=bar(1:days,[rates(1:days,1) rates(1:days,2)]);
        qqq(1).FaceColor='flat'; qqq(2).FaceColor='flat';
        qqq(1).CData=[hitColorPlot];qqq(2).CData=[faColorPlot];
        ylim([0 1]);title(['H and FA rate ']); xlabel('days');
%         subplot(1,3,3);www=bar(1:2,[mean(rates(lastdaynum-3:lastdaynum,1)) mean(rates(lastdaynum-3:lastdaynum,2))]);
%         www.FaceColor='flat'; www.CData(1:2,:)=[hitcolor;facolor];
%         title(['H and FA rate last 3 days']); xticklabels({'hit','fa'});
        
        if savefig
            cd(pathsave);
            saveas(perfFig,[subjlist{nbsubj} '-LearningSummary-' num2str(nbins) '.' filetype]);
            saveas(perfFig,[subjlist{nbsubj} '-LearningSummary-' num2str(nbins) '.png']);
            close(perfFig);
        end
    else
    end
    
%% TO PLOT OPTO
% if nbsubj==1
%     optoplot=1;
% elseif nbsubj==6
%     optoplot=1; 
% elseif nbsubj==7
%     optoplot=1; 
% else
%     optoplot=0;
% end
optoplot=1;

% bar graphs for opto
if optoplot==1 % now make bar graphs, averaged, for all conditions 
    
% rates(i,:) = [rhit rfa phit pfa ohit ofa ...
%          rfhit rffa cofffa conffa phit2 pfa2 ...
%          othit otfa ochit ocfa];
% rates 1 -2 reinforced H and FA
% rates 3-4 probe
% rates 5-6 opto full trial
% rates 7 - 14 stuff we dont need for now (catch light off, catch light on)
% rates 15 - 16 - tone
% rates 17 - 18 - choice 

    eeeFig=figure('Position', [10 10 725 575]);hold on;

%     if contains(path,'1')==1        
%         if nbsubj ==4 % cohort 1
%             mgbDays = [1 1 0 0 1 1 0 1 1 1];mgbDays=logical(mgbDays);
%             icDays = [0 0 1 1 0 0 1 0 0 0];icDays=logical(icDays);
%             expRange=1:8; % Cohort 1
%         elseif nbsubj==5 % sk163
%             mgbDays = [1 1 0 0 1 1 0 0 1 1 1 0 0];mgbDays=logical(mgbDays);
%             icDays = [0 0 1 1 0 0 1 1 0 0 0 1 1];icDays=logical(icDays);
%             expRange=1:8; % Cohort 1
%         else
%             mgbDays = [1 1 0 0 1 1 0 0 1 1 1];mgbDays=logical(mgbDays);
%             icDays = [0 0 1 1 0 0 1 1 0 0 0]; icDays=logical(icDays);
%             expRange=1:8; % Cohort 1
%         end
%     elseif contains(path,'2')==1
%         if nbsubj==1
%             expRange=25:32; % cohort 2
%             mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
%                 1 1 0 0 1 1 0 0 ...
%                 0 0 0 ];mgbDays=logical(mgbDays);
%             icDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 1 1 0 0 1 1 ...
%                 0 0 0];icDays=logical(icDays);
%         elseif nbsubj==2
%             expRange=22:32; % cohort 2
%             mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 0 1 1 0 0 1 1 0 0 ...
%                 0 0 0];mgbDays=logical(mgbDays);
%             icDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
%                 1 1 1 0 0 1 1 0 0 1 1 ...
%                 0 0 0];icDays=logical(icDays);
%         end
%     elseif contains(path,'3')==1
%         if nbsubj==1 %sk176
%             mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 0 0 0 0 0 0 0 ...
%                 1 0 1 1 0 1 1 0 1 0 0 1];
%             mgbDays=logical(mgbDays);
%             icDays=[0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 0 0 0 0 0 0 0 ... 
%                 0 1 0 0 1 0 0 1 0 1 1 0];
%             icDays=logical(icDays);
%             expRange=20:length(mgbDays); % cohort 3 sk177
%         elseif nbsubj==2 %sk177
%             mgbDays=[0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 0];
%             mgbDays=logical(mgbDays);
%             icDays=[0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1 1];
%             icDays=logical(icDays);
%             expRange=8:17; % cohort 3 sk177
%         elseif nbsubj==3 %sk178
%             mgbDays= [0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 0 0 0 0 0 0 0];
%             mgbDays=logical(mgbDays);
%             icDays=[0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 0 0 0 0 0 0 0];
%             icDays=logical(icDays);
%             expRange=15:17;
%         elseif nbsubj==4 %sk179
%             mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 0 0 1 1 1 0 0 1 0 0];
%             mgbDays=logical(mgbDays);
%             icDays = [0 0 0 0 0 0 0 0 0 0 ...
%                 0 0 0 0 0 0 0 1 1 0 1 1];
%             icDays=logical(icDays);
%             expRange=15:22;  % cohort 3 sk179
%         end
    if contains(path,'4')==1
        if nbsubj==2 %sk183
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 1 1 0 0 1 1 0 0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 1 1 0 0 1 1];
            icDays=logical(icDays);
            expRange=13:length(mgbDays); % cohort 3 sk183
        elseif nbsubj==5
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 1];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0];
            icDays=logical(icDays);
            expRange=length(mgbDays); % cohort 3 sk186   
        elseif nbsubj==8
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 1 1 0 0 1 1 1 0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 1 1 0 0 0 1];
            icDays=logical(icDays);
            expRange=13:length(mgbDays); % cohort 3 sk189        
        end
    elseif contains(path,'5')==1
        if nbsubj==1 %sk190 
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ... %d30
                0 0 1 1 0 0 0 1 1 0 0 0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 1 ...
                1 0 0 0 0 1 1 0 0 0 0 0]; % as of 10/3, need to update
            icDays=logical(icDays);
            expRange=21:length(mgbDays);
        elseif nbsubj==2 %sk191
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 1 1 0 0 0 ...
                0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 1 1 0 0 0 0 1 1 0 ...
                0]; % mgb, ic catch, lob
            icDays=logical(icDays);
            expRange=22:length(mgbDays);
        elseif nbsubj==5 % 194 good
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 0 1 1 0 ...
                0 1 0 0 0 0 0]; % catch ic, catch ic, catch mgb, lob
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 1 1 0 0 1 1 0 0 1 ...
                1 0 1 1 0 0 0]; % catch ic, catch ic, catch mgb, lob
            icDays=logical(icDays);
            expRange=12:length(icDays);
        elseif nbsubj==6 % 195
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 1 0 0 0 ...
                0 1 0 0]; % until 10/3
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 1 1 1 ...
                1 0 1 1]; 
            icDays=logical(icDays);
            expRange=24:length(icDays);  
        elseif nbsubj==7 %sk196 good 
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 1 0 0 1 1 ... 
                1 0]; % missing day 9/24 aka d26
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 1 1 0 1 1 0 0 ...
                0 0];
            icDays=logical(icDays);
            expRange=23:length(mgbDays);        
        else
            mgbDays = [0 0 0 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 0 0 0]; icDays=logical(icDays);
            expRange=1:5; 
        end
        
    elseif contains(path,'6')==1
        if nbsubj==1 % sk198
            mgbDays=[0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 1 1 ...
                0 1 0 1 0 0 0 ...
                0 0 0]; %... % catch mgb, catch ic, lob
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                1 0 1 0 1 0 1 ...
                0 0 0]; %...
            icDays=logical(icDays);
            expRange=38:length(mgbDays);
        elseif nbsubj==6 %sk203, started on d9
            mgbDays=[0 0 0 0 0 0 0 0 ...
                1 1 0 0 0 1 1 0 0];%...
%                 1 0 0]; % catch MGB, catch IC, LOB
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 0 1 1];%...
%                 0 1 0]; % catch MGB, catch IC, LOB
            icDays=logical(icDays);
            expRange=9:length(mgbDays);
        elseif nbsubj==7 %sk204, started on d29
            mgbDays=[0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ... %15
                0 0 0 0 0 0 0 0 0 ...
                0 0 ... %26
                0 0 1 1 0 1 0 1 ...
                0 0 ...  % no opto, then ic tone
                0 0 0 0]; % catch ic, catch mgb, lob
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ...
                0 0 ...
                1 0 0 0 1 0 1 0 ...
                0 1 ...
                0 0 1 0]; % catch ic, catch mgb, ic tone, lob
            icDays=logical(icDays);
            expRange=27:length(mgbDays); 
        else
            mgbDays = [0 0 0 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 0 0 0]; icDays=logical(icDays);
            expRange=1:5; 
        end
    elseif contains(path,'IC_cohort_1')==1
        if nbsubj==1 % sk218 d30 start
            mgbDays=[0 0 0 0 0 0 ...
                0 0 0 0 0 0 ...
                0 0 0 0 0 0 ...
                0 0 0 0 0 0 ...
                0 0 ...
                0 0 0 0 0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 ...
                0 0 0 0 0 0 ...
                0 0 0 0 0 0 ...
                0 0 0 0 0 0 ...
                0 1 ...
                1 1 1 1 1];
            icDays=logical(icDays);
            expRange=28:length(icDays);
        elseif nbsubj==2 %sk219
            mgbDays=[0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 0];
            icDays=logical(icDays);
            expRange=6:length(icDays);
        elseif nbsubj==3 %sk220
            mgbDays=[0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 0];
            icDays=logical(icDays);
            expRange=7:length(icDays);
    
        end
        
    end
    
    mgbDays=mgbDays(expRange);
    icDays=icDays(expRange);

    reinfcolor= [0.4,0.4,0.4];
    optocolor=[102/255 178/255 255/255];

%     % make a grid of two by three, where each row is a targeted area (MGB or
    % IC) and each column is a type of opto (full trial/tone/choice)
    if sum(mgbDays)==0
        subplot(3, 1, 1);eee=bar([mean(rates(expRange(icDays),1)) nanmean(rates(expRange(icDays),5)); ...
            mean(rates(expRange(icDays),2)) nanmean(rates(expRange(icDays),6))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor]; hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),5)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),6)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); ylabel('Average Action Rate'); title(['IC Full Opto']);
        % legend('light off','light on');

        subplot(3,1,2);eee=bar([mean(rates(expRange(icDays),1)) nanmean(rates(expRange(icDays),15)); ...
            mean(rates(expRange(icDays),2)) nanmean(rates(expRange(icDays),16))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),15)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),16)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('IC Tone Opto');
        % legend('light off','light on');

        subplot(3,1,3);eee=bar([mean(rates(expRange(icDays),1)) nanmean(rates(expRange(icDays),17)); ...
            mean(rates(expRange(icDays),2)) nanmean(rates(expRange(icDays),18))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),17)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),18)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('IC Choice Opto');
        % legend('light off','light on');

        if savefig
            cd(pathsave);
            saveas(eeeFig,[subjlist{nbsubj} '-OptoBarScatter-' num2str(nbins) '.' filetype]);
            saveas(eeeFig,[subjlist{nbsubj} '-OptoBarScatter-' num2str(nbins) '.png']);
            close(eeeFig);
        end
        
        optomeanMat{nbsubj+nbsubj,1}=subjlist(nbsubj);
        optomeanMat{nbsubj+nbsubj,2}=mean(rates(:,1));
        optomeanMat{nbsubj+nbsubj,3}=mean(rates(:,2)); 
        optomeanMat{nbsubj+nbsubj,4}=nanmean(rates(:,5));
        optomeanMat{nbsubj+nbsubj,5}=nanmean(rates(:,6));
        optomeanMat{nbsubj+nbsubj,6}=nanmean(rates(:,15));
        optomeanMat{nbsubj+nbsubj,7}=nanmean(rates(:,16));
        optomeanMat{nbsubj+nbsubj,8}=mean(rates(:,17));
        optomeanMat{nbsubj+nbsubj,9}=mean(rates(:,18));
        optomeanMat{nbsubj+nbsubj,10}=mean(rates(expRange(mgbDays),1)); % MGB Data
        optomeanMat{nbsubj+nbsubj,11}=mean(rates(expRange(mgbDays),2));
        optomeanMat{nbsubj+nbsubj,12}=nanmean(rates(expRange(mgbDays),5));
        optomeanMat{nbsubj+nbsubj,13}=nanmean(rates(expRange(mgbDays),6));
        optomeanMat{nbsubj+nbsubj,14}=nanmean(rates(expRange(mgbDays),15));
        optomeanMat{nbsubj+nbsubj,15}=nanmean(rates(expRange(mgbDays),16));
        optomeanMat{nbsubj+nbsubj,16}=mean(rates(expRange(mgbDays),17));
        optomeanMat{nbsubj+nbsubj,17}=mean(rates(expRange(mgbDays),18));
        optomeanMat{nbsubj+nbsubj,18}=mean(rates(expRange(icDays),1)); % IC Data
        optomeanMat{nbsubj+nbsubj,19}=mean(rates(expRange(icDays),2));
        optomeanMat{nbsubj+nbsubj,20}=nanmean(rates(expRange(icDays),5));
        optomeanMat{nbsubj+nbsubj,21}=nanmean(rates(expRange(icDays),6));
        optomeanMat{nbsubj+nbsubj,22}=nanmean(rates(expRange(icDays),15));
        optomeanMat{nbsubj+nbsubj,23}=nanmean(rates(expRange(icDays),16));
        optomeanMat{nbsubj+nbsubj,24}=mean(rates(expRange(icDays),17));
        optomeanMat{nbsubj+nbsubj,25}=mean(rates(expRange(icDays),18));

        optomeanMat{nbsubj+nbsubj+1,1}=subjlist(nbsubj);
        optomeanMat{nbsubj+nbsubj+1,2}=rates(:,1);
        optomeanMat{nbsubj+nbsubj+1,3}=rates(:,2); 
        optomeanMat{nbsubj+nbsubj+1,4}=rates(:,5);
        optomeanMat{nbsubj+nbsubj+1,5}=rates(:,6);
        optomeanMat{nbsubj+nbsubj+1,6}=rates(:,15);
        optomeanMat{nbsubj+nbsubj+1,7}=rates(:,16);
        optomeanMat{nbsubj+nbsubj+1,8}=rates(:,17);
        optomeanMat{nbsubj+nbsubj+1,9}=rates(:,18);
        optomeanMat{nbsubj+nbsubj+1,10}=rates(expRange(mgbDays),1); % MGB Data
        optomeanMat{nbsubj+nbsubj+1,11}=rates(expRange(mgbDays),2);
        optomeanMat{nbsubj+nbsubj+1,12}=rates(expRange(mgbDays),5);
        optomeanMat{nbsubj+nbsubj+1,13}=rates(expRange(mgbDays),6);
        optomeanMat{nbsubj+nbsubj+1,14}=rates(expRange(mgbDays),15);
        optomeanMat{nbsubj+nbsubj+1,15}=rates(expRange(mgbDays),16);
        optomeanMat{nbsubj+nbsubj+1,16}=rates(expRange(mgbDays),17);
        optomeanMat{nbsubj+nbsubj+1,17}=rates(expRange(mgbDays),18);
        optomeanMat{nbsubj+nbsubj+1,18}=rates(expRange(icDays),1); % IC Data
        optomeanMat{nbsubj+nbsubj+1,19}=rates(expRange(icDays),2);
        optomeanMat{nbsubj+nbsubj+1,20}=rates(expRange(icDays),5);
        optomeanMat{nbsubj+nbsubj+1,21}=rates(expRange(icDays),6);
        optomeanMat{nbsubj+nbsubj+1,22}=rates(expRange(icDays),15);
        optomeanMat{nbsubj+nbsubj+1,23}=rates(expRange(icDays),16);
        optomeanMat{nbsubj+nbsubj+1,24}=rates(expRange(icDays),17);
        optomeanMat{nbsubj+nbsubj+1,25}=rates(expRange(icDays),18);
        optomeanMat{nbsubj+nbsubj+1,26}=MAT{nbsubj,1};
        optomeanMat{nbsubj+nbsubj+1,27}=rates;
        
    else
        subplot(3,3,1);eee=bar([mean(rates(expRange,1)) nanmean(rates(expRange,5)); ...
            mean(rates(expRange,2)) nanmean(rates(expRange,6))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange,1),1),1),(rates(expRange,1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange,1),1),2),(rates(expRange,2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange,1),1),1),(rates(expRange,5)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange,1),1),2),(rates(expRange,6)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); ylabel('Average Action Rate'); legend('light off','light on');
        title([subjlist{nbsubj} ' (' expnames{explist(nbsubj)} ') ' 'Full Trial Opto']);

        subplot(3,3,2);eee=bar([mean(rates(expRange,1)) nanmean(rates(expRange,15)); ...
            mean(rates(expRange,2)) nanmean(rates(expRange,16))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange,1),1),1),(rates(expRange,1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange,1),1),2),(rates(expRange,2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange,1),1),1),(rates(expRange,15)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange,1),1),2),(rates(expRange,16)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('Tone Only Opto');
    %     legend('light off','light on');

        subplot(3,3,3);eee=bar([mean(rates(expRange,1)) nanmean(rates(expRange,17)); ...
            mean(rates(expRange,2)) nanmean(rates(expRange,18))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor]; hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange,1),1),1),(rates(expRange,1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange,1),1),2),(rates(expRange,2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange,1),1),1),(rates(expRange,17)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange,1),1),2),(rates(expRange,18)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('Choice Only Opto');
    %     legend('light off','light on');
    % 

        subplot(3,3,4);eee=bar([mean(rates(expRange(mgbDays),1)) nanmean(rates(expRange(mgbDays),5)); ...
            mean(rates(expRange(mgbDays),2)) nanmean(rates(expRange(mgbDays),6))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(mgbDays),1),1),1),(rates(expRange(mgbDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(mgbDays),1),1),2),(rates(expRange(mgbDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(mgbDays),1),1),1),(rates(expRange(mgbDays),5)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(mgbDays),1),1),2),(rates(expRange(mgbDays),6)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); ylabel('Average Action Rate'); title(['MGB Full Opto']);
        % legend('light off','light on');

        subplot(3,3,5);eee=bar([mean(rates(expRange(mgbDays),1)) nanmean(rates(expRange(mgbDays),15)); ...
            mean(rates(expRange(mgbDays),2)) nanmean(rates(expRange(mgbDays),16))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor]; hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(mgbDays),1),1),1),(rates(expRange(mgbDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(mgbDays),1),1),2),(rates(expRange(mgbDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(mgbDays),1),1),1),(rates(expRange(mgbDays),15)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(mgbDays),1),1),2),(rates(expRange(mgbDays),16)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('MGB Tone Opto');
        % legend('light off','light on');

        subplot(3,3,6);eee=bar([mean(rates(expRange(mgbDays),1)) nanmean(rates(expRange(mgbDays),17)); ...
            mean(rates(expRange(mgbDays),2)) nanmean(rates(expRange(mgbDays),18))]); %hit, full opto hit, fa, opto fa 
         eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(mgbDays),1),1),1),(rates(expRange(mgbDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(mgbDays),1),1),2),(rates(expRange(mgbDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(mgbDays),1),1),1),(rates(expRange(mgbDays),17)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(mgbDays),1),1),2),(rates(expRange(mgbDays),18)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('MGB Choice Opto');
        % legend('light off','light on');

        subplot(3,3,7);eee=bar([mean(rates(expRange(icDays),1)) nanmean(rates(expRange(icDays),5)); ...
            mean(rates(expRange(icDays),2)) nanmean(rates(expRange(icDays),6))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor]; hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),5)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),6)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); ylabel('Average Action Rate'); title(['IC Full Opto']);
        % legend('light off','light on');


        subplot(3,3,8);eee=bar([mean(rates(expRange(icDays),1)) nanmean(rates(expRange(icDays),15)); ...
            mean(rates(expRange(icDays),2)) nanmean(rates(expRange(icDays),16))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),15)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),16)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('IC Tone Opto');
        % legend('light off','light on');


        subplot(3,3,9);eee=bar([mean(rates(expRange(icDays),1)) nanmean(rates(expRange(icDays),17)); ...
            mean(rates(expRange(icDays),2)) nanmean(rates(expRange(icDays),18))]); %hit, full opto hit, fa, opto fa
        eee(1).FaceColor='flat'; eee(2).FaceColor='flat'; eee(1).CData=[reinfcolor;reinfcolor];hold on;
        scatter(repmat(eee(1).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),1)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(eee(1).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),2)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        eee(2).CData=[optocolor;optocolor];ylim([0 1]);
        scatter(repmat(eee(2).XEndPoints(1),size(rates(expRange(icDays),1),1),1),(rates(expRange(icDays),17)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(eee(2).XEndPoints(2),size(rates(expRange(icDays),1),1),2),(rates(expRange(icDays),18)),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        xticklabels({'hit','fa'}); title('IC Choice Opto');
        % legend('light off','light on');

        if savefig
            cd(pathsave);
            saveas(eeeFig,[subjlist{nbsubj} '-OptoBarScatter-' num2str(nbins) '.' filetype]);
            saveas(eeeFig,[subjlist{nbsubj} '-OptoBarScatter-' num2str(nbins) '.png']);
            close(eeeFig);
        end

        optomeanMat{nbsubj+nbsubj,1}=subjlist(nbsubj);
        optomeanMat{nbsubj+nbsubj,2}=mean(rates(:,1));
        optomeanMat{nbsubj+nbsubj,3}=mean(rates(:,2)); 
        optomeanMat{nbsubj+nbsubj,4}=nanmean(rates(:,5));
        optomeanMat{nbsubj+nbsubj,5}=nanmean(rates(:,6));
        optomeanMat{nbsubj+nbsubj,6}=nanmean(rates(:,15));
        optomeanMat{nbsubj+nbsubj,7}=nanmean(rates(:,16));
        optomeanMat{nbsubj+nbsubj,8}=mean(rates(:,17));
        optomeanMat{nbsubj+nbsubj,9}=mean(rates(:,18));
        optomeanMat{nbsubj+nbsubj,10}=mean(rates(expRange(mgbDays),1)); % MGB Data
        optomeanMat{nbsubj+nbsubj,11}=mean(rates(expRange(mgbDays),2));
        optomeanMat{nbsubj+nbsubj,12}=nanmean(rates(expRange(mgbDays),5));
        optomeanMat{nbsubj+nbsubj,13}=nanmean(rates(expRange(mgbDays),6));
        optomeanMat{nbsubj+nbsubj,14}=nanmean(rates(expRange(mgbDays),15));
        optomeanMat{nbsubj+nbsubj,15}=nanmean(rates(expRange(mgbDays),16));
        optomeanMat{nbsubj+nbsubj,16}=mean(rates(expRange(mgbDays),17));
        optomeanMat{nbsubj+nbsubj,17}=mean(rates(expRange(mgbDays),18));
        optomeanMat{nbsubj+nbsubj,18}=mean(rates(expRange(icDays),1)); % IC Data
        optomeanMat{nbsubj+nbsubj,19}=mean(rates(expRange(icDays),2));
        optomeanMat{nbsubj+nbsubj,20}=nanmean(rates(expRange(icDays),5));
        optomeanMat{nbsubj+nbsubj,21}=nanmean(rates(expRange(icDays),6));
        optomeanMat{nbsubj+nbsubj,22}=nanmean(rates(expRange(icDays),15));
        optomeanMat{nbsubj+nbsubj,23}=nanmean(rates(expRange(icDays),16));
        optomeanMat{nbsubj+nbsubj,24}=mean(rates(expRange(icDays),17));
        optomeanMat{nbsubj+nbsubj,25}=mean(rates(expRange(icDays),18));

        optomeanMat{nbsubj+nbsubj+1,1}=subjlist(nbsubj);
        optomeanMat{nbsubj+nbsubj+1,2}=rates(:,1);
        optomeanMat{nbsubj+nbsubj+1,3}=rates(:,2); 
        optomeanMat{nbsubj+nbsubj+1,4}=rates(:,5);
        optomeanMat{nbsubj+nbsubj+1,5}=rates(:,6);
        optomeanMat{nbsubj+nbsubj+1,6}=rates(:,15);
        optomeanMat{nbsubj+nbsubj+1,7}=rates(:,16);
        optomeanMat{nbsubj+nbsubj+1,8}=rates(:,17);
        optomeanMat{nbsubj+nbsubj+1,9}=rates(:,18);
        optomeanMat{nbsubj+nbsubj+1,10}=rates(expRange(mgbDays),1); % MGB Data
        optomeanMat{nbsubj+nbsubj+1,11}=rates(expRange(mgbDays),2);
        optomeanMat{nbsubj+nbsubj+1,12}=rates(expRange(mgbDays),5);
        optomeanMat{nbsubj+nbsubj+1,13}=rates(expRange(mgbDays),6);
        optomeanMat{nbsubj+nbsubj+1,14}=rates(expRange(mgbDays),15);
        optomeanMat{nbsubj+nbsubj+1,15}=rates(expRange(mgbDays),16);
        optomeanMat{nbsubj+nbsubj+1,16}=rates(expRange(mgbDays),17);
        optomeanMat{nbsubj+nbsubj+1,17}=rates(expRange(mgbDays),18);
        optomeanMat{nbsubj+nbsubj+1,18}=rates(expRange(icDays),1); % IC Data
        optomeanMat{nbsubj+nbsubj+1,19}=rates(expRange(icDays),2);
        optomeanMat{nbsubj+nbsubj+1,20}=rates(expRange(icDays),5);
        optomeanMat{nbsubj+nbsubj+1,21}=rates(expRange(icDays),6);
        optomeanMat{nbsubj+nbsubj+1,22}=rates(expRange(icDays),15);
        optomeanMat{nbsubj+nbsubj+1,23}=rates(expRange(icDays),16);
        optomeanMat{nbsubj+nbsubj+1,24}=rates(expRange(icDays),17);
        optomeanMat{nbsubj+nbsubj+1,25}=rates(expRange(icDays),18);
        optomeanMat{nbsubj+nbsubj+1,26}=MAT{nbsubj,1};
        optomeanMat{nbsubj+nbsubj+1,27}=rates;
    end
else
end

save('summaryData.mat','optomeanMat');
disp('saved opto data to mat file.');
end
    
%% single animal, single day bar plot to see hit/m/cr/fa 

% figure;hold on;bbb=bar([sum(reinf&hit) sum(reinf)/2-sum(reinf&hit) sum(reinf)/2-sum(reinf&fa) sum(reinf&fa)...
%     sum(probe&hit) sum(probe)/2-sum(probe&hit) sum(probe)/2-sum(probe&fa) sum(probe&fa)...
%     sum(opto&hit) sum(opto)/2-sum(opto&hit) sum(opto)/2-sum(opto&fa) sum(opto&fa)...
%     sum(copto&hit) sum(copto)/2-sum(copto&hit) sum(copto)/2-sum(copto&fa) sum(copto&fa)],'FaceColor','flat');
% xticks([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);
% title([subjlist{nbsubj} ' (' expnames{explist(nbsubj)} ') - Raw Data']);
% xticklabels({'Reinf Hit','Reinf Miss','Reinf CR','Reinf FA','Probe Hit','Probe Miss','Probe CR','Probe FA',...
%     'Full Opto Hit','Full Opto miss','Full Opto CR','Full Opto FA',...
%     'Choice Opto Hit','Choice Opto Miss','Choice Opto CR','Choice Opto FA'});
% bbb.CData(1:4,:)=[0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5;0.5 0.5 0.5]; % Reinf colors
% bbb.CData(5:8,:)=[.7 .7 .7;.7 .7 .7;.7 .7 .7;.7 .7 .7]; % probe colors
% bbb.CData(9:12,:)=[0,0,255/255;0,0,255/255;0,0,255/255;0,0,255/255]; % Opto colors
% bbb.CData(13:16,:)=[102/255 178/255 255/255;102/255 178/255 255/255;102/255 178/255 255/255;102/255 178/255 255/255]; % Choice pto colors

   %% now make plots averaged across test and ctl for MGB
   ctl = find(explist==2);
%    ctl=ctl(2:4);
   test = find(explist==1);
   ggg=figure;
   subplot(3,2,1);ggg=bar([mean(cell2mat(optomeanMat(ctl+1,9))) nanmean(cell2mat(optomeanMat(ctl+1,11)))...
       mean(cell2mat(optomeanMat(ctl+1,10))) nanmean(cell2mat(optomeanMat(ctl+1,12)))]); %hit, full opto hit, fa, opto fa
    ggg(1).FaceColor='flat';  ggg(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ggg(1).XEndPoints(1),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,9))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(2),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,11))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ggg(1).XEndPoints(3),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,10))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(4),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,12))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('CTL MGB Full Trial ');
    legend('','light off','','light on');
    
    subplot(3,2,2);ggg=bar([mean(cell2mat(optomeanMat(test+1,9))) nanmean(cell2mat(optomeanMat(test+1,11)))...
       mean(cell2mat(optomeanMat(test+1,10))) nanmean(cell2mat(optomeanMat(test+1,12)))]); %hit, full opto hit, fa, opto fa
    ggg(1).FaceColor='flat';  ggg(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ggg(1).XEndPoints(1),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,9))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(2),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,11))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ggg(1).XEndPoints(3),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,10))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(4),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,12))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('Test MGB Full Trial ');  
    
    subplot(3,2,3);ggg=bar([mean(cell2mat(optomeanMat(ctl+1,9))) nanmean(cell2mat(optomeanMat(ctl+1,15)))...
       mean(cell2mat(optomeanMat(ctl+1,10))) nanmean(cell2mat(optomeanMat(ctl+1,16)))]); %hit, full opto hit, fa, opto fa
    ggg(1).FaceColor='flat';  ggg(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ggg(1).XEndPoints(1),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,9))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEnsubdPoints(2),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,15))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ggg(1).XEndPoints(3),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,10))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(4),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,16))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('CTL MGB Choice');
    
    subplot(3,2,4);ggg=bar([mean(cell2mat(optomeanMat(test+1,9))) nanmean(cell2mat(optomeanMat(test+1,15)))...
       mean(cell2mat(optomeanMat(test+1,10))) nanmean(cell2mat(optomeanMat(ctl+1,16)))]); %hit, full opto hit, fa, opto fa
    ggg(1).FaceColor='flat';  ggg(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ggg(1).XEndPoints(1),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,9))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(2),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,15))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ggg(1).XEndPoints(3),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,10))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(4),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,16))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''});  title('Test MGB Choice ');      
    
   subplot(3,2,5);ggg=bar([mean(cell2mat(optomeanMat(ctl+1,9))) nanmean(cell2mat(optomeanMat(ctl+1,13)))...
       mean(cell2mat(optomeanMat(ctl+1,10))) nanmean(cell2mat(optomeanMat(ctl+1,14)))]); %hit, full opto hit, fa, opto fa
    ggg(1).FaceColor='flat';  ggg(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ggg(1).XEndPoints(1),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,9))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(2),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,13))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ggg(1).XEndPoints(3),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,10))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(4),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,14))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('CTL MGB Tone');
    
    subplot(3,2,6);ggg=bar([mean(cell2mat(optomeanMat(test+1,9))) nanmean(cell2mat(optomeanMat(test+1,13)))...
       mean(cell2mat(optomeanMat(test+1,10))) nanmean(cell2mat(optomeanMat(ctl+1,14)))]); %hit, full opto hit, fa, opto fa
    ggg(1).FaceColor='flat';  ggg(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ggg(1).XEndPoints(1),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,9))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(2),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,13))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ggg(1).XEndPoints(3),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,10))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ggg(1).XEndPoints(4),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,14))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''});  title('Test MGB Tone ');          
   
 %% now make the same plot but for the IC
 
   ooo=figure;
   subplot(3,2,1);ooo=bar([mean(cell2mat(optomeanMat(ctl+1,17))) nanmean(cell2mat(optomeanMat(ctl+1,19)))...
       mean(cell2mat(optomeanMat(ctl+1,18))) nanmean(cell2mat(optomeanMat(ctl+1,20)))]); %hit, full opto hit, fa, opto fa
    ooo(1).FaceColor='flat';  ooo(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ooo(1).XEndPoints(1),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,17))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(2),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,19))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ooo(1).XEndPoints(3),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,18))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(4),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,20))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('CTL IC Full Trial ');
    legend('','light off','','light on');
    
    subplot(3,2,2);ooo=bar([mean(cell2mat(optomeanMat(test+1,17))) nanmean(cell2mat(optomeanMat(test+1,19)))...
       mean(cell2mat(optomeanMat(test+1,18))) nanmean(cell2mat(optomeanMat(test+1,20)))]); %hit, full opto hit, fa, opto fa
    ooo(1).FaceColor='flat';  ooo(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ooo(1).XEndPoints(1),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,17))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(2),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,19))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ooo(1).XEndPoints(3),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,18))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(4),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,20))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('Test IC Full Trial ');  
    
    subplot(3,2,3);ooo=bar([mean(cell2mat(optomeanMat(ctl+1,17))) nanmean(cell2mat(optomeanMat(ctl+1,23)))...
       mean(cell2mat(optomeanMat(ctl+1,18))) nanmean(cell2mat(optomeanMat(ctl+1,24)))]); %hit, full opto hit, fa, opto fa
    ooo(1).FaceColor='flat';  ooo(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ooo(1).XEndPoints(1),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,17))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(2),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,23))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ooo(1).XEndPoints(3),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,18))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(4),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,24))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('CTL IC Choice');
    
    subplot(3,2,4);ooo=bar([mean(cell2mat(optomeanMat(test+1,17))) nanmean(cell2mat(optomeanMat(test+1,23)))...
       mean(cell2mat(optomeanMat(test+1,18))) nanmean(cell2mat(optomeanMat(test+1,24)))]); %hit, full opto hit, fa, opto fa
    ooo(1).FaceColor='flat';  ooo(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ooo(1).XEndPoints(1),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,17))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(2),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,23))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ooo(1).XEndPoints(3),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,18))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(4),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,24))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''});  title('Test IC Choice ');      
    
   subplot(3,2,5);ooo=bar([mean(cell2mat(optomeanMat(ctl+1,17))) nanmean(cell2mat(optomeanMat(ctl+1,21)))...
       mean(cell2mat(optomeanMat(ctl+1,18))) nanmean(cell2mat(optomeanMat(ctl+1,22)))]); %hit, full opto hit, fa, opto fa
    ooo(1).FaceColor='flat';  ooo(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ooo(1).XEndPoints(1),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,17))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(2),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,21))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ooo(1).XEndPoints(3),size(cell2mat(optomeanMat(ctl+1,1)),1),1),(cell2mat(optomeanMat(ctl+1,18))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(4),size(cell2mat(optomeanMat(ctl+1,1)),1),2),(cell2mat(optomeanMat(ctl+1,22))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''}); title('CTL IC Tone');
    
    subplot(3,2,6);ooo=bar([mean(cell2mat(optomeanMat(test+1,17))) nanmean(cell2mat(optomeanMat(test+1,21)))...
       mean(cell2mat(optomeanMat(test+1,18))) nanmean(cell2mat(optomeanMat(test+1,22)))]); %hit, full opto hit, fa, opto fa
    ooo(1).FaceColor='flat';  ooo(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];hold on;
    scatter(repmat(ooo(1).XEndPoints(1),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,17))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(2),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,21))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    scatter(repmat(ooo(1).XEndPoints(3),size(cell2mat(optomeanMat(test+1,1)),1),1),(cell2mat(optomeanMat(test+1,18))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    scatter(repmat(ooo(1).XEndPoints(4),size(cell2mat(optomeanMat(test+1,1)),1),2),(cell2mat(optomeanMat(test+1,22))),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
    xticklabels({'        hit','','      fa',''});  title('Test IC Tone ');          
   

 %% Plot catch trials
% 
% % explist = [2 1 2 1 1 2 1 2]; % 1 = test, 2 = ctl   subjlist = {'CD077','CD078','CD079','CD080','CD081','CD082','CD083','CD084'};
% 
% % explist = [2 1 2 2 1 1 2 1 2 1 1 2]; % 1 = test, 2 = ctl
% cond = 'all';
% % explist = [2 1 2 2 1 1 2 1 0 1 0 2]; % only no learner removed
% % cond = 'nonlearnerremoved';
% % explist = [2 1 2 2 1 0 2 1 0 1 0 0]; % all bad removed
% % cond = 'badremoved';
% % cd060 (test): good then bad
% % cd065 (test): no learning at all
% % cd063 (ctl): no learning, only one day with perf>1
% % cd066 (ctl): good then bad
ctl = find(explist==2);
test = find(explist==1);
savefig = false;

fig=figure; % FA rate
c = 1; t = 9;
for i=1:nSubj    
    if ismember(i,ctl)
        subplot(2,8,c);hold on;
        nFiles = length(subjrates{i}(:,9));
        plot(1:nFiles,subjrates{i}(:,9),'k.-','markersize',10);
        plot(1:nFiles,subjrates{i}(:,10),'b.-','markersize',10);
        c=c+1;
        ylim([0 1]);ylabel('FA rate');
        title(subjlist{i});
        subplot(2,4,4);hold on;
        plot(1:nFiles,subjrates{i}(:,9),'k.-','markersize',10);
        plot(1:nFiles,subjrates{i}(:,10),'b.-','markersize',10);
        ylim([0 1]);
    elseif ismember(i,test)
        subplot(2,8,t);hold on;
        nFiles = length(subjrates{i}(:,9));
        plot(1:nFiles,subjrates{i}(:,9),'k.-','markersize',10);
        plot(1:nFiles,subjrates{i}(:,10),'b.-','markersize',10);
        t=t+1;
        ylim([0 1]);ylabel('FA rate');  
        title(subjlist{i});
        subplot(2,4,8);hold on;
        plot(1:nFiles,subjrates{i}(:,9),'k.-','markersize',10);
        plot(1:nFiles,subjrates{i}(:,10),'b.-','markersize',10);
        ylim([0 1]);
    end
end
subplot(2,8,7);hold on; % plot catch with light on the brain
data = [];groups = [];
for i=1:nSubj
    if ismember(i,ctl)
        data = [data;subjrates{i}(1:end-1,9);subjrates{i}(1:end-1,10)];
        groups = [groups;ones(nFiles-1,1);2*ones(nFiles-1,1)];
    end
end
anovabar(data,groups);hold on;
plot([1 2],[data(groups==1) data(groups==2)], 'k.-');
[~,p] = ttest(data(groups==1),data(groups==2));
title(['Ttest paired, p=' num2str(p)]);
ylim([0 1]);
subplot(2,8,8);hold on; % plot catch with light outside the brain
data = [];groups = [];
for i=1:nSubj
    if ismember(i,ctl)
        data = [data;subjrates{i}(end,9);subjrates{i}(end,10)];
        groups = [groups;ones(1,1);2*ones(1,1)];
    end
end
anovabar(data,groups);hold on;
plot(groups, data, 'k.-');
[~,p] = ttest(data(groups==1),data(groups==2));
title(['Ttest paired, p=' num2str(p)]);
ylim([0 1]);

subplot(2,8,15);hold on;
data = [];groups = [];
for i=1:nSubj
    if ismember(i,test)
        data = [data;subjrates{i}(1:end-1,9);subjrates{i}(1:end-1,10)];
        groups = [groups;ones(nFiles-1,1);2*ones(nFiles-1,1)];
    end
end
% anovabar(data,groups);hold on;
% [~,p] = ttest(data(groups==1),data(groups==2));
% title(['Ttest paired, p=' num2str(p)]);
% plot([1 2],[data(groups==1) data(groups==2)], 'b.-');
% ylim([0 1]);
% subplot(2,8,16);hold on;
% data = [];groups = [];
% for i=1:nSubj
%     if ismember(i,test)
%         data = [data;subjrates{i}(end,9);subjrates{i}(end,10)];
%         groups = [groups;ones(1,1);2*ones(1,1)];
%     end
% end
% anovabar(data,groups);hold on;
% [~,p] = ttest(data(groups==1),data(groups==2));
% title(['Ttest paired, p=' num2str(p)]);
% plot(groups, data, 'b.-');
% ylim([0 1]);
% 
% if savefig
%     cd(pathsave);
%     saveas(fig,['CatchTrials-' cond '.pdf']);
%     close(fig);
% end
% 
% %% Plot Lick histograms during catch trials
% % explist = [2 1 2 2 1 1 2 1 2 1 1 2]; % 1 = test, 2 = ctl
% % cond = 'all';
% % explist = [2 1 2 2 1 1 2 1 0 1 0 2]; % only no learner removed
% % cond = 'nonlearnerremoved';
% % explist = [2 1 2 2 1 0 2 1 0 1 0 0]; % all bad removed
% % cond = 'badremoved';
% % cd060 (test): good then bad
% % cd065 (test): no learning at all
% % cd063 (ctl): no learning, only one day with perf>1
% % cd066 (ctl): good then bad
% ctl = find(explist==2);
% test = find(explist==1);
% savefig = false;
% 
% fig=figure; % lick hist (relative: start - stim end - dead - resp wind end)
% c = 1; t = 7;
% x = linspace(-1,4,nbins);
% smoothfactor = 6;
% ylimm = [0 3];
% for i=1:nSubj    
%     if ismember(i,ctl)
%         subplot(2,6,c);hold on;
%         plot(x,Smooth(sum(lickhistcOFF{i}),smoothfactor),'k-');
%         plot(x,Smooth(sum(lickhistcON{i}),smoothfactor),'b-');
%         c=c+1;
%         ylim(ylimm);
%         PlotIntervals([0 0.1]);
%         title(subjlist{i});
%         subplot(2,6,6);hold on;
%         plot(x,Smooth(sum(lickhistcOFF{i}),smoothfactor),'k-');
%         plot(x,Smooth(sum(lickhistcON{i}),smoothfactor),'b-');
%         ylim(ylimm);
%         PlotIntervals([0 0.1]);
%     elseif ismember(i,test)
%         subplot(2,6,t);hold on;
%         plot(x,Smooth(sum(lickhistcOFF{i}),smoothfactor),'k-');
%         plot(x,Smooth(sum(lickhistcON{i}),smoothfactor),'b-');
%         t=t+1;
%         ylim(ylimm);
%         PlotIntervals([0 0.1]);
%         title(subjlist{i});
%         subplot(2,6,12);hold on;
%         plot(x,Smooth(sum(lickhistcOFF{i}),smoothfactor),'k-');
%         plot(x,Smooth(sum(lickhistcON{i}),smoothfactor),'b-');
%         ylim(ylimm);
%         PlotIntervals([0 0.1]);
%     end
% end
% 
% if savefig
%     cd(pathsave);
%     saveas(fig,['Average-lickPSTH-CatchTrials-' cond '.pdf']);
%     close(fig);
% end
% 

%% Plot difference performance light-off light-on

ctl = find(explist==2);
test = find(explist==1);
maxday = max(matrix(:,SESS));
days = 1:maxday;
savefig = false;

fig=figure;hold on;
subplot(1,3,1);hold on;
diffrodprime_test = odprime(test,days)-rdprime(test,days);
% diffrodprime_ctl = odprime(test,days)-rdprime(test,days);
plot(1:maxday,diffrodprime_test','-o','color',[0.5 0.7 0.9],'LineWidth',2);
% plot(1:maxday,diffrodprime_ctl','k-');
PlotHVLines(0,'h','k');
shadedErrorBar(1:maxday,nanmean(diffrodprime_test),nansem(diffrodprime_test),{'color','b','linewidth',2},0.5);
% shadedErrorBar(1:maxday,nanmean(diffrodprime_ctl),nansem(diffrodprime_ctl),{'color','k','linewidth',2},0.5);
xlabel('Days');ylabel("d' reinforced light on - light off");
xlim([1 maxday]);
subplot(1,3,2);hold on; 
plot(rdprime(test,days),diffrodprime_test,'.','color',[0.5 0.7 0.9],'markersize',20);
% plot(rdprime(ctl,days),diffrodprime_ctl,'k.','markersize',10);
xlabel("d' reinforced light-off");
ylabel("d' reinforced light on - light off");
plot(-1:5,0*ones(length(-1:5)),'k');
plot(0*ones(length(-1:5)),-5:1,'k');
ylim([-5 1]);xlim([-1 5]);
% title(['ctl n=' num2str(length(ctl)) ', test n=' num2str(length(test))]);

subplot(1,3,3);hold on; 
plot(rdprime(test,days),odprime(test,days),'.','color',[0.5 0.7 0.9],'markersize',20);
% plot(rdprime(ctl,days),odprime(ctl,days),'k.','markersize',10);
xlabel("d' reinforced light-off");
ylabel("d' reinforced light-on");
PlotHVLines(0,'v','k');
PlotHVLines(0,'h','k');
ylim([0 5]);xlim([0 5]);
plot(0:5,0:5,'k');

rdprime_test = rdprime(test,days);rdprime_test = rdprime_test(:);
odprime_test = odprime(test,days);odprime_test = odprime_test(:);
nop = any(isnan([rdprime_test odprime_test]),2);
[rho_test,yfit] = LinearRegression(rdprime_test(~nop),odprime_test(~nop));
plot(rdprime_test(~nop),yfit,'b--');
[~,p_test] = corr(rdprime_test(~nop),odprime_test(~nop)); 
% rdprime_ctl = rdprime(ctl,days);rdprime_ctl = rdprime_ctl(:);
% odprime_ctl = odprime(ctl,days);odprime_ctl = odprime_ctl(:);
% nop = any(isnan([rdprime_ctl odprime_ctl]),2);
% [rho_ctl,yfit] = LinearRegression(rdprime_ctl(~nop),odprime_ctl(~nop)); 
% plot(rdprime_ctl(~nop),yfit,'k--');
% [~,p_ctl] = corr(rdprime_ctl(~nop),odprime_ctl(~nop)); 
% title(['r2 ctl = ' num2str(rho_ctl) ', test = ' num2str(rho_test)]);

if savefig
    cd(pathsave);
    saveas(fig,['PerformanceComparison-lightOnOff-' cond '.pdf']);
    close(fig);
end

%%


ctl = find(explist==1);
test = find(explist==1);
maxday = 17;
days = 1:maxday;
savefig = false;

fig=figure;hold on;
subplot(2,2,1);hold on;
diffrodprime_test = odprime(test,days)-rdprime(test,days);
diffrodprime_ctl = odprime(ctl,days)-rdprime(ctl,days);
plot(1:maxday,diffrodprime_test','b-');
plot(1:maxday,diffrodprime_ctl','k-');
PlotHVLines(0,'h','k');
shadedErrorBar(1:maxday,nanmean(diffrodprime_test),nansem(diffrodprime_test),{'color','b','linewidth',2},0.5);
shadedErrorBar(1:maxday,nanmean(diffrodprime_ctl),nansem(diffrodprime_ctl),{'color','k','linewidth',2},0.5);
xlabel('Days');ylabel("d' reinforced light on - light off");
xlim([1 maxday]);
subplot(2,2,2);hold on; 
plot(rdprime(test,days),diffrodprime_test,'b.','markersize',10);
plot(rdprime(ctl,days),diffrodprime_ctl,'k.','markersize',10);
xlabel("d' reinforced light-off");
ylabel("d' reinforced light on - light off");
plot(-1:5,0*ones(length(-1:5)),'k');
plot(0*ones(length(-1:5)),-5:1,'k');
ylim([-5 1]);xlim([-1 5]);
title(['ctl n=' num2str(length(ctl)) ', test n=' num2str(length(test))]);

subplot(2,2,3);hold on; 
plot(rdprime(test,days),odprime(test,days),'b.','markersize',10);
plot(rdprime(ctl,days),odprime(ctl,days),'k.','markersize',10);
xlabel("d' reinforced light-off");
ylabel("d' reinforced light-on");
PlotHVLines(0,'v','k');
PlotHVLines(0,'h','k');
ylim([0 5]);xlim([0 5]);
plot(0:5,0:5,'k');

rdprime_test = rdprime(test,days);rdprime_test = rdprime_test(:);
odprime_test = odprime(test,days);odprime_test = odprime_test(:);
nop = any(isnan([rdprime_test odprime_test]),2);
[rho_test,yfit] = LinearRegression(rdprime_test(~nop),odprime_test(~nop));
plot(rdprime_test(~nop),yfit,'b--');
[~,p_test] = corr(rdprime_test(~nop),odprime_test(~nop)); 
rdprime_ctl = rdprime(ctl,days);rdprime_ctl = rdprime_ctl(:);
odprime_ctl = odprime(ctl,days);odprime_ctl = odprime_ctl(:);
nop = any(isnan([rdprime_ctl odprime_ctl]),2);
[rho_ctl,yfit] = LinearRegression(rdprime_ctl(~nop),odprime_ctl(~nop)); 
plot(rdprime_ctl(~nop),yfit,'k--');
[~,p_ctl] = corr(rdprime_ctl(~nop),odprime_ctl(~nop)); 
title(['r2 ctl = ' num2str(rho_ctl) ', test = ' num2str(rho_test)]);

if savefig
    cd(pathsave);
    saveas(fig,['PerformanceComparison-lightOnOff-' cond '.pdf']);
    close(fig);
end

%% Comparison light-on vs light-off end of learning

% Pool variables to analyze

maxdays = 20; % used to be 30 
rhit = nan(nSubj,maxdays);
rfa = nan(nSubj,maxdays);
phit = nan(nSubj,maxdays);
pfa = nan(nSubj,maxdays);
ohit = nan(nSubj,maxdays);
ofa = nan(nSubj,maxdays);
cofffa = nan(nSubj,maxdays);
confa = nan(nSubj,maxdays);
% 
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
% 
% 
% Compute percent correct
rpc = (rhit+(1-rfa))/2*100; 
opc = (ohit+(1-ofa))/2*100; 
% 
% starts = find(sum(isnan(opc))~=nSubj);
% starts(end)=[]; % remove light outside the brain
% 
% colors = {'c','y','b','g'};
% 
% figure;
% for s=1:nSubj
%     subplot(2,2,s);hold on;
%     plot([rpc(s,starts)' opc(s,starts)']','k.-');
%     xlim([0.5 2.5]);
%     set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
%     ylim([0 100]);ylabel('Proportion correct');
%     [~,p] = ttest(rpc(s,starts)',opc(s,starts)');
%     title([subjlist{s} ', p=' num2str(p)]);
% end
% 
% figure;hold on;
% subplot(2,2,1);hold on;
% lightoff = rpc(:,starts)';
% lighton = opc(:,starts)';
% plot([lightoff(:) lighton(:)]','.-','color',[0.8 0.8 0.8]);
% plot([1 2],[nanmean(lightoff(:)) nanmean(lighton(:))],'k.','markersize',15);
% xlim([0.5 2.5]);
% set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
% ylim([0 100]);ylabel('Proportion correct');
% subplot(2,2,2);hold on;
% for s=1:nSubj
%     plot([rpc(s,starts)' opc(s,starts)']','-','color',colors{s});
% end
% xlim([0.5 2.5]);
% set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
% ylim([0 100]);ylabel('Proportion correct');
% plot([1 2],[nanmean(lightoff(:)) nanmean(lighton(:))],'k.','markersize',15);
% if sum([jbtest(lightoff(:));jbtest(lighton(:))]) == 0
%     [~,p] = ttest(lightoff(:),lighton(:));
%     title(['paired t-test, p=' num2str(p)]);
% else
%     p = signrank(lightoff(:),lighton(:));
%     title(['paired Wilcoxon, p=' num2str(p)]);
% end
% 
% subplot(2,2,3);hold on;
% plot([1 2],[nanmean(lightoff(:)) nanmean(lighton(:))],'k.','markersize',15);
% for s=1:nSubj
%     lightoff = nanmean(rpc(s,starts),2);
%     lighton = nanmean(opc(s,starts),2);
%     plot([lightoff lighton]','.-','color',colors{s});
% end
% xlim([0.5 2.5]);
% set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
% ylim([0 100]);ylabel('Proportion correct');
% 
% subplot(2,2,4);hold on;
% title('Light outside the brain');
% outside = find(sum(isnan(opc))~=nSubj);
% outside = outside(end);
% lightoff = rpc(:,outside)';
% lighton = opc(:,outside)';
% plot([lightoff(:) lighton(:)]','k.-');
% plot([1 2],[nanmean(lightoff(:)) nanmean(lighton(:))],'b.','markersize',15);
% xlim([0.5 2.5]);
% set(gca,'xtick',[1 2],'xticklabel',{'light-off','light-on'});
% ylim([0 100]);ylabel('Proportion correct');

%% Group plots
% subjlist = {'CD055','CD056','CD057','CD058','CD059','CD060','CD061','CD062','CD063','CD064','CD065','CD066'};
% explist = [2 1 2 2 1 1 2 1 2 1 1 2]; % 1 = test, 2 = ctl
% explist = [2 1 2 2 1 1 2 1 2 1 1 2 1 1 1 2]; % 1 = test, 2 = ctl

ctl = find(explist==2);
test = find(explist==1);

reinfTestVsCtl = true;
optoTestVsCtl = true;
probeTestVsCtl = true;
reinfVsOptoTest = true;
reinfVsOptoCtl = true;
maxFiles = 15; %used to be set to 30 

if reinfTestVsCtl
    
    % Plot REIFORCED data, TEST vs CTL 
    fig = figure; hold on;
    subplot(2,4,1); hold on; % HIT rate, test vs control
    title('REINFORCED, TEST vs CTL');
    colors = {'k','b'};
    nFiles = maxFiles;
    ctl_subjrates_hit = nan(nSubj/2,nFiles);
    ctl_subjrates_fa = nan(nSubj/2,nFiles);
    test_subjrates_hit = nan(nSubj/2,nFiles);
    test_subjrates_fa = nan(nSubj/2,nFiles);
    c = 1 ; t = 1;
    for i=1:nSubj    
        if ismember(i,ctl)
            nFiles = length(subjrates{i}(:,1)');
            ctl_subjrates_hit(c,1:nFiles) = subjrates{i}(:,1)';
            ctl_subjrates_fa(c,1:nFiles) = subjrates{i}(:,2)';
            c=c+1;
        elseif ismember(i,test)
            nFiles = length(subjrates{i}(:,1)');
            test_subjrates_hit(t,1:nFiles) = subjrates{i}(:,1)';
            test_subjrates_fa(t,1:nFiles) = subjrates{i}(:,2)';
            t=t+1;
        end       
    end
    nFiles = maxFiles;
    plot(1:nFiles,test_subjrates_hit','.-','color',colors{2});
    plot(1:nFiles,ctl_subjrates_hit','.-','color',colors{1});
    ylabel('HIT rate');xlabel('Days');
    ctl_m = nanmean(ctl_subjrates_hit);
    ctl_s = nansem(ctl_subjrates_hit);
%     shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'color','k'});
    test_m = nanmean(test_subjrates_hit);
    test_s = nansem(test_subjrates_hit);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'-','color','b'});

    subplot(2,4,2); hold on; % FA rate, test vs control
    plot(1:nFiles,test_subjrates_fa','.--','color',colors{2});
    plot(1:nFiles,ctl_subjrates_fa','.--','color',colors{1});
    test_m = nanmean(test_subjrates_fa);
    test_s = nansem(test_subjrates_fa);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'color','k'});
    ctl_m = nanmean(ctl_subjrates_fa);
    ctl_s = nansem(ctl_subjrates_fa);
%     shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'-','color','b'});
    ylabel('FA rate');xlabel('Days');

    subplot(2,4,3); hold on; % correct
    c=1; t=1;
    nFiles = maxFiles;
    correct_ctl = nan(nSubj/2,nFiles);
    correct_test = nan(nSubj/2,nFiles);
    for i=1:nSubj
        correct = (subjrates{i}(:,1) + (1-subjrates{i}(:,2)) ) /2; %.*subjnbctxtoutcome{i};
        if ismember(i,ctl)
            nFiles = length(correct');
            correct_ctl(c,1:nFiles) = correct';
            c=c+1;
        elseif ismember(i,test)
            nFiles = length(correct');
            correct_test(t,1:nFiles) = correct';
            t=t+1;
        end
    end
    nFiles = maxFiles;
    plot(1:nFiles,correct_test','.-','color',colors{2});
    plot(1:nFiles,correct_ctl','.-','color',colors{1});
    ylim([0 1]);
    ylabel('Portion correct'); xlabel('Days');

    subplot(2,4,4); hold on; % d primes
    plot(1:nFiles,rdprime(ctl,:),'k.-');
    plot(1:nFiles,rdprime(test,:),'b.-');
%     shadedErrorBar(1:nFiles,mean(rdprime(ctl,:)),sem(rdprime(ctl,:)),{'color','k'});
%     shadedErrorBar(1:nFiles,mean(rdprime(test,:)),sem(rdprime(test,:)),{'color','b'});
    ylim([-1 5]);
    ylabel("d'"); xlabel('Days');
    
    subplot(2,4,5); hold on; % HIT lick latency, test vs control
    plot(1:size(mr_hit(ctl,:)',1),mr_hit(ctl,:)','k.-');
    plot(1:size(mr_hit(test,:)',1),mr_hit(test,:)','b.-');
    ylabel('HIT lick latency (s)');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,4,6); hold on; % FA lick latency, test vs control
    plot(1:size(mr_fa(ctl,:)',1),mr_fa(ctl,:)','k.--');
    plot(1:size(mr_fa(test,:)',1),mr_fa(test,:)','b.--');
    ylabel('FA lick latency (s)');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,4,7); hold on; % HIT lick latency, test vs control
    plot(1:size(mlr_hit(ctl,:)',1),mlr_hit(ctl,:)','k.-');
    plot(1:size(mlr_hit(test,:)',1),mlr_hit(test,:)','b.-');
    ylabel('HIT lick rate (Hz)');xlabel('Days');
    ylim([0 10]);

    subplot(2,4,8); hold on; % FA lick latency, test vs control
    plot(1:size(mlr_fa(ctl,:)',1),mlr_fa(ctl,:)','k.--');
    plot(1:size(mlr_fa(test,:)',1),mlr_fa(test,:)','b.--');
    ylabel('FA lick rate (Hz)');xlabel('Days');
    ylim([0 10]);
    
    fig = figure; hold on;
    subplot(2,4,1); hold on; % HIT rate, test vs control
    title('REINFORCED, TEST vs CTL');
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    data = [test_subjrates_hit(:);ctl_subjrates_hit(:)];
    barwitherrn(data,groups);hold on;
    g1 = linspace(0.75,1.25,nSubj/2*nFiles);
    g2 = linspace(1.75,2.25,nSubj/2*nFiles);
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel('HIT rate');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,4,2); hold on; % FA rate, test vs control
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    data = [test_subjrates_fa(:);ctl_subjrates_fa(:)];
    barwitherrn(data,groups);hold on;
    g1 = linspace(0.75,1.25,nSubj/2*nFiles);
    g2 = linspace(1.75,2.25,nSubj/2*nFiles);
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel('FA rate');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,4,3); hold on; % correct
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    data = [correct_test(:);correct_ctl(:)];
    barwitherrn(data,groups);hold on;
    g1 = linspace(0.75,1.25,nSubj/2*nFiles);
    g2 = linspace(1.75,2.25,nSubj/2*nFiles);
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel('Portion correct');xlabel('Days');
    ylim([0 1]);
    
    subplot(2,4,4); hold on; % d primes
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    ctl_mr_hit = rdprime(ctl,:);
    test_mr_hit = rdprime(test,:);
    data = [test_mr_hit(:);ctl_mr_hit(:)];
    anovabar(data,groups);hold on;
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel("d'");xlabel('Days');
    ylim([-1 4]);
    h1 = jbtest(data(groups==1));h2 = jbtest(data(groups==2));
    if h1+h2 ==0
        disp('Data are normal');
    else
        p = signrank(data(groups==1),data(groups==2));
        title(['p=' num2str(p)]);
    end

    subplot(2,4,5); hold on; % HIT lick latency, test vs control
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    ctl_mr_hit = mr_hit(ctl,:);
    test_mr_hit = mr_hit(test,:);
    data = [test_mr_hit(:);ctl_mr_hit(:)];
    barwitherrn(data,groups);hold on;
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel('HIT lick latency (s)');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,4,6); hold on; % FA lick latency, test vs control
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    ctl_mr_fa = mr_fa(ctl,:);
    test_mr_fa = mr_fa(test,:);
    data = [test_mr_fa(:);ctl_mr_fa(:)];
    barwitherrn(data,groups);hold on;
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel('FA lick latency (s)');xlabel('Days');
    ylim([0 1.5]);
    
    subplot(2,4,7); hold on; % HIT lick rate, test vs control
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    ctl_mr_hit = mlr_hit(ctl,:);
    test_mr_hit = mlr_hit(test,:);
    data = [test_mr_hit(:);ctl_mr_hit(:)];
    barwitherrn(data,groups);hold on;
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel('HIT lick rate (Hz)');xlabel('Days');
    ylim([0 10]);

    subplot(2,4,8); hold on; % FA lick rate, test vs control
    groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
    ctl_mr_fa = mlr_fa(ctl,:);
    test_mr_fa = mlr_fa(test,:);
    data = [test_mr_fa(:);ctl_mr_fa(:)];
    barwitherrn(data,groups);hold on;
    plot(g1,data(groups==1),'o');
    plot(g2,data(groups==2),'o');
    set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
    ylabel('FA lick rate (Hz)');xlabel('Days');
    ylim([0 10]);

    
end
    
if optoTestVsCtl

    % Plot OPTO data, TEST vs CTL 
    ctl = find(explist==2);
    test = find(explist==1);

    fig = figure; hold on;
    subplot(2,3,1); hold on; % HIT rate, test vs control
    title('OPTO trials, TEST vs CTL');
    colors = {'k','b'};
    ctl_subjrates_hit = nan(nSubj/2,nFiles);
    ctl_subjrates_fa = nan(nSubj/2,nFiles);
    test_subjrates_hit = nan(nSubj/2,nFiles);
    test_subjrates_fa = nan(nSubj/2,nFiles);
    c = 1 ; t = 1;
    for i=1:nSubj    
        if ismember(i,ctl)
            nFiles = length(subjrates{i}(:,5)');
            ctl_subjrates_hit(c,1:nFiles) = subjrates{i}(:,5)';
            ctl_subjrates_fa(c,1:nFiles) = subjrates{i}(:,6)';
            c=c+1;
        elseif ismember(i,test)
            nFiles = length(subjrates{i}(:,5)');
            test_subjrates_hit(t,1:nFiles) = subjrates{i}(:,5)';
            test_subjrates_fa(t,1:nFiles) = subjrates{i}(:,6)';
            t=t+1;
        end       
    end
    nFiles = maxFiles;    
    plot(1:nFiles,test_subjrates_hit','.-','color',colors{2});
    plot(1:nFiles,ctl_subjrates_hit','.-','color',colors{1});
    ylabel('HIT rate');xlabel('Days');
    ctl_m = nanmean(ctl_subjrates_hit);
    ctl_s = nansem(ctl_subjrates_hit);
%     shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'color','k'});
    test_m = nanmean(test_subjrates_hit);
    test_s = nansem(test_subjrates_hit);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'-','color','b'});

    subplot(2,3,2); hold on; % FA rate, test vs control
    plot(1:nFiles,test_subjrates_fa','.--','color',colors{2});
    plot(1:nFiles,ctl_subjrates_fa','.--','color',colors{1});
    test_m = nanmean(test_subjrates_fa);
    test_s = nansem(test_subjrates_fa);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'color','k'});
    ctl_m = nanmean(ctl_subjrates_fa);
    ctl_s = nansem(ctl_subjrates_fa);
%     shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'-','color','b'});
    ylabel('FA rate');xlabel('Days');

    subplot(2,3,3); hold on; % correct
    c=1; t=1;
    correct_ctl = nan(nSubj/2,nFiles);
    correct_test = nan(nSubj/2,nFiles);
    for i=1:nSubj
        correct = (subjrates{i}(:,5) + (1-subjrates{i}(:,6)) ) /2; %.*subjnbctxtoutcome{i};
        if ismember(i,ctl)
            nFiles = length(correct');
            correct_ctl(c,1:nFiles) = correct';
            c=c+1;
        elseif ismember(i,test)
            nFiles = length(correct');
            correct_test(t,1:nFiles) = correct';
            t=t+1;
        end
    end
    nFiles = maxFiles;
    plot(1:nFiles,correct_test','.-','color',colors{2});
    plot(1:nFiles,correct_ctl','.-','color',colors{1});
    ylim([0 1]);
    ylabel('Portion correct'); xlabel('Days');

    nFiles = size(mo_hit(ctl,:),2);
    subplot(2,3,4); hold on; % HIT lick latency, test vs control
    plot(1:nFiles,mo_hit(ctl,:)','k.-');
    plot(1:nFiles,mo_hit(test,:)','b.-');
    ylabel('HIT lick latency (s)');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,3,5); hold on; % FA lick latency, test vs control
    plot(1:nFiles,mo_fa(ctl,:)','k.--');
    plot(1:nFiles,mo_fa(test,:)','b.--');
    ylabel('FA lick latency (s)');xlabel('Days');
    ylim([0 1.5]);

    nFiles = maxFiles;
    subplot(2,3,6); hold on; % d primes
    plot(1:nFiles,odprime(ctl,:)','.-','color',colors{1});
    plot(1:nFiles,odprime(test,:)','.-','color',colors{2});
    ylim([-1 5]);
    ylabel("d'"); xlabel('Days');

end

if probeTestVsCtl

    % Plot PROBE data, TEST vs CTL 
    ctl = find(explist==2);
    test = find(explist==1);

    fig = figure; hold on;
    subplot(2,3,1); hold on; % HIT rate, test vs control
    title('PROBE, TEST vs CTL');
    colors = {'k','b'};
    ctl_subjrates_hit = nan(nSubj/2,nFiles);
    ctl_subjrates_fa = nan(nSubj/2,nFiles);
    test_subjrates_hit = nan(nSubj/2,nFiles);
    test_subjrates_fa = nan(nSubj/2,nFiles);
    c = 1 ; t = 1;
    for i=1:nSubj    
        if ismember(i,ctl)
            nFiles = length(subjrates{i}(:,3)');
            ctl_subjrates_hit(c,1:nFiles) = subjrates{i}(:,3)';
            ctl_subjrates_fa(c,1:nFiles) = subjrates{i}(:,4)';
            c=c+1;
        elseif ismember(i,test)
            nFiles = length(subjrates{i}(:,3)');
            test_subjrates_hit(t,1:nFiles) = subjrates{i}(:,3)';
            test_subjrates_fa(t,1:nFiles) = subjrates{i}(:,4)';
            t=t+1;
        end       
    end
    nFiles = maxFiles;
%     plot(1:nFiles,test_subjrates_hit','.-','color',colors{2});
%     plot(1:nFiles,ctl_subjrates_hit','.-','color',colors{1});
    ylabel('HIT rate');xlabel('Days');
    ctl_m = nanmean(ctl_subjrates_hit);
    ctl_s = nansem(ctl_subjrates_hit);
    shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'color','k'});
    test_m = nanmean(test_subjrates_hit);
    test_s = nansem(test_subjrates_hit);
    shadedErrorBar(1:nFiles,test_m,test_s,{'-','color','b'});

    subplot(2,3,2); hold on; % FA rate, test vs control
%     plot(1:nFiles,test_subjrates_fa','.--','color',colors{2});
%     plot(1:nFiles,ctl_subjrates_fa','.--','color',colors{1});
    test_m = nanmean(test_subjrates_fa);
    test_s = nansem(test_subjrates_fa);
    shadedErrorBar(1:nFiles,test_m,test_s,{'color','k'});
    ctl_m = nanmean(ctl_subjrates_fa);
    ctl_s = nansem(ctl_subjrates_fa);
    shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'-','color','b'});
    ylabel('FA rate');xlabel('Days');

    subplot(2,3,3); hold on; % correct
    c=1; t=1;
    correct_ctl = nan(nSubj/2,nFiles);
    correct_test = nan(nSubj/2,nFiles);
    for i=1:nSubj
        correct = (subjrates{i}(:,3) + (1-subjrates{i}(:,4)) ) /2; %.*subjnbctxtoutcome{i};
        if ismember(i,ctl)
            nFiles = length(correct');
            correct_ctl(c,1:nFiles) = correct';
            c=c+1;
        elseif ismember(i,test)
            nFiles = length(correct');
            correct_test(t,1:nFiles) = correct';
            t=t+1;
        end
    end
    nFiles = maxFiles;
    plot(1:nFiles,correct_test','.-','color',colors{2});
    plot(1:nFiles,correct_ctl','.-','color',colors{1});
    ylim([0 1]);
    ylabel('Portion correct'); xlabel('Days');

    subplot(2,3,4); hold on; % HIT lick latency, test vs control
    plot(1:size(mp_hit(ctl,:)',1),mp_hit(ctl,:)','k.-');
    plot(1:size(mp_hit(test,:)',1),mp_hit(test,:)','b.-');
    ylabel('HIT lick latency (s)');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,3,5); hold on; % FA lick latency, test vs control
    plot(1:size(mp_fa(test,:)',1),mp_fa(ctl,:)','k.--');
    plot(1:size(mp_fa(test,:)',1),mp_fa(test,:)','b.--');
    ylabel('FA lick latency (s)');xlabel('Days');
    ylim([0 1.5]);

    subplot(2,3,6); hold on; % d primes
%     plot(1:nFiles,pdprime(ctl,:)','.-','color',colors{1});
%     plot(1:nFiles,pdprime(test,:)','.-','color',colors{2});
    shadedErrorBar(1:nFiles,mean(pdprime(ctl,:)),sem(pdprime(ctl,:)),{'color','k'});
    shadedErrorBar(1:nFiles,mean(pdprime(test,:)),sem(pdprime(test,:)),{'color','b'});
    ylim([-1 5]);
    ylabel("d'"); xlabel('Days');

end

if reinfVsOptoTest

    % Plot REINF vs OPTO data, TEST
    test = find(explist==1);

    fig = figure; hold on;
    subplot(2,3,1); hold on; % HIT rate, test vs control
    title('REINF vs OPTO trials, TEST');
    colors = {'k','b'};
    test_subjrates_hit_r = nan(nSubj/2,nFiles);
    test_subjrates_fa_r = nan(nSubj/2,nFiles);
    test_subjrates_hit_o = nan(nSubj/2,nFiles);
    test_subjrates_fa_o = nan(nSubj/2,nFiles);
    c = 1 ; t = 1;
    for i=1:nSubj    
        if ismember(i,test)
            nFiles = length(subjrates{i}(:,1)');
            test_subjrates_hit_r(c,1:nFiles) = subjrates{i}(:,1)';
            test_subjrates_fa_r(c,1:nFiles) = subjrates{i}(:,2)';
            c=c+1;
            test_subjrates_hit_o(t,1:nFiles) = subjrates{i}(:,5)';
            test_subjrates_fa_o(t,1:nFiles) = subjrates{i}(:,6)';
            t=t+1;
        end       
    end
    nFiles = maxFiles;
    plot(1:nFiles,test_subjrates_hit_r','.-','color',colors{1});
    plot(1:nFiles,test_subjrates_hit_o','.-','color',colors{2});
    ylabel('HIT rate');xlabel('Days');    
    test_m = nanmean(test_subjrates_hit_r);
    test_s = nansem(test_subjrates_hit_r);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'color','k'});
    test_m = nanmean(test_subjrates_hit_o);
    test_s = nansem(test_subjrates_hit_o);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'-','color','b'});
    ylim([0 1]);
    
    subplot(2,3,2); hold on; % FA rate 
    plot(1:nFiles,test_subjrates_fa_r','.--','color',colors{1});
    plot(1:nFiles,test_subjrates_fa_o','.--','color',colors{2});
    test_m = nanmean(test_subjrates_fa_r);
    test_s = nansem(test_subjrates_fa_r);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'color','k'});
    ctl_m = nanmean(test_subjrates_fa_o);
    ctl_s = nansem(test_subjrates_fa_o);
%     shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'-','color','b'});
    ylabel('FA rate');xlabel('Days');
    ylim([0 1]);
    
    subplot(2,3,3); hold on; % correct
    t=1;
    correct_test_r = nan(nSubj/2,nFiles);
    correct_test_o = nan(nSubj/2,nFiles);
    for i=1:nSubj
        if ismember(i,test)
            nFiles = length(((subjrates{i}(:,1) + (1-subjrates{i}(:,2)) ) /2 )');
            correct_test_r(t,1:nFiles) = ((subjrates{i}(:,1) + (1-subjrates{i}(:,2)) ) /2 )';
            correct_test_o(t,1:nFiles) = ((subjrates{i}(:,5) + (1-subjrates{i}(:,6)) ) /2 )';
            t=t+1;
        end
    end
    nFiles = maxFiles;
    plot(1:nFiles,correct_test_o','.-','color',colors{2});
    plot(1:nFiles,correct_test_r','.-','color',colors{1});
    ylim([0 1]);
    ylabel('Portion correct'); xlabel('Days');

    subplot(2,3,4); hold on; % HIT lick latency 
    plot(1:size(mr_hit(test,:)',1),mr_hit(test,:)','k.-');
    plot(1:size(mo_hit(test,:)',1),mo_hit(test,:)','b.-');
    ylabel('HIT lick latency (s)');xlabel('Days');
    ylim([0 2]);

    subplot(2,3,5); hold on; % FA lick latency 
    plot(1:size(mr_fa(test,:)',1),mr_fa(test,:)','k.--');
    plot(1:size(mo_fa(test,:)',1),mo_fa(test,:)','b.--');
    ylabel('FA lick latency (s)');xlabel('Days');
    ylim([0 2]);

    subplot(2,3,6); hold on; % d primes
    plot(1:nFiles,rdprime(test,:)','.-','color',colors{1});
    plot(1:nFiles,odprime(test,:)','.-','color',colors{2});
    ylim([-1 5]);
    ylabel("d'"); xlabel('Days');

end

if reinfVsOptoCtl

    % Plot REINF vs OPTO data, CTL
    test = find(explist==1);

    fig = figure; hold on;
    subplot(2,3,1); hold on; % HIT rate, test vs control
    title('REINF vs OPTO trials, CTL');
    colors = {'k','b'};
    ctl_subjrates_hit_r = nan(nSubj/2,nFiles);
    ctl_subjrates_fa_r = nan(nSubj/2,nFiles);
    ctl_subjrates_hit_o = nan(nSubj/2,nFiles);
    ctl_subjrates_fa_o = nan(nSubj/2,nFiles);
    c = 1 ; t = 1;
    for i=1:nSubj    
        if ismember(i,ctl)
            nFiles = length(subjrates{i}(:,1)');
            ctl_subjrates_hit_r(c,1:nFiles) = subjrates{i}(:,1)';
            ctl_subjrates_fa_r(c,1:nFiles) = subjrates{i}(:,2)';
            c=c+1;
            ctl_subjrates_hit_o(t,1:nFiles) = subjrates{i}(:,5)';
            ctl_subjrates_fa_o(t,1:nFiles) = subjrates{i}(:,6)';
            t=t+1;
        end       
    end
    nFiles = maxFiles;
    ylabel('HIT rate');xlabel('Days');
    plot(1:nFiles,ctl_subjrates_hit_r','.-','color',colors{1});
    plot(1:nFiles,ctl_subjrates_hit_o','.-','color',colors{2});
    ctl_m = nanmean(ctl_subjrates_hit_r);
    ctl_s = nansem(ctl_subjrates_hit_r);
%     shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'color','k'});
    test_m = nanmean(ctl_subjrates_hit_o);
    test_s = nansem(ctl_subjrates_hit_o);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'-','color','b'});

    subplot(2,3,2); hold on; % FA rate 
    plot(1:nFiles,ctl_subjrates_fa_r','.--','color',colors{1});
    plot(1:nFiles,ctl_subjrates_fa_o','.--','color',colors{2});
    test_m = nanmean(ctl_subjrates_fa_r);
    test_s = nansem(ctl_subjrates_fa_r);
%     shadedErrorBar(1:nFiles,test_m,test_s,{'color','k'});
    ctl_m = nanmean(ctl_subjrates_fa_o);
    ctl_s = nansem(ctl_subjrates_fa_o);
%     shadedErrorBar(1:nFiles,ctl_m,ctl_s,{'-','color','b'});
    ylabel('FA rate');xlabel('Days');

    subplot(2,3,3); hold on; % correct
    c=1;
    correct_ctl_r = nan(nSubj/2,nFiles);
    correct_ctl_o = nan(nSubj/2,nFiles);
    for i=1:nSubj
        if ismember(i,ctl)
            nFiles = length(((subjrates{i}(:,1) + (1-subjrates{i}(:,2)) ) /2 )');
            correct_ctl_r(c,1:nFiles) = ((subjrates{i}(:,1) + (1-subjrates{i}(:,2)) ) /2 )';
            correct_ctl_o(c,1:nFiles) = ((subjrates{i}(:,5) + (1-subjrates{i}(:,6)) ) /2 )';
            c=c+1;
        end
    end
    nFiles = maxFiles;
    plot(1:nFiles,correct_ctl_o','.-','color',colors{2});
    plot(1:nFiles,correct_ctl_r','.-','color',colors{1});
    ylim([0 1]);
    ylabel('Portion correct'); xlabel('Days');

    subplot(2,3,4); hold on; % HIT lick latency 
    plot(1:size(mr_hit(ctl,:)',1),mr_hit(ctl,:)','k.-');
    plot(1:size(mo_hit(ctl,:)',1),mo_hit(ctl,:)','b.-');
    ylabel('HIT lick latency (s)');xlabel('Days');
    ylim([0 2]);

    subplot(2,3,5); hold on; % FA lick latency 
    plot(1:size(mr_fa(ctl,:)',1),mr_fa(ctl,:)','k.--');
    plot(1:size(mo_fa(ctl,:)',1),mo_fa(ctl,:)','b.--');
    ylabel('FA lick latency (s)');xlabel('Days');
    ylim([0 2]);

    subplot(2,3,6); hold on; % d primes
    plot(1:nFiles,rdprime(ctl,:)','.-','color',colors{1});
    plot(1:nFiles,odprime(ctl,:)','.-','color',colors{2});
    ylim([-1 5]);
    ylabel("d'"); xlabel('Days');

%     fig = figure; hold on;
%     subplot(2,3,1); hold on; % HIT rate, test vs control
%     title('REINF vs OPTO trials, CTL');
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [ctl_subjrates_hit_r(:);ctl_subjrates_hit_o(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = ctl_subjrates_hit_r(i,:);
%         g2 = ctl_subjrates_hit_o(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     ylabel('HIT rate');xlabel('Days');
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylim([0 1.5]);
% 
%     subplot(2,3,2); hold on; % FA rate 
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [ctl_subjrates_fa_r(:);ctl_subjrates_fa_o(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = ctl_subjrates_fa_r(i,:);
%         g2 = ctl_subjrates_fa_o(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     ylabel('FA rate');xlabel('Days');
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylabel('FA rate');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,3); hold on; % correct
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [correct_ctl_r(:);correct_ctl_o(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = correct_ctl_r(i,:);
%         g2 = correct_ctl_o(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylabel('Portion correct');xlabel('Days');
%     ylim([0 1]);
% 
%     subplot(2,3,4); hold on; % HIT lick latency 
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mr_hit = mr_hit(ctl,:);
%     ctl_mo_hit = mo_hit(ctl,:);
%     data = [ctl_mr_hit(:);ctl_mo_hit(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = ctl_mr_hit(i,:);
%         g2 = ctl_mo_hit(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     ylabel('HIT lick latency (s)');xlabel('Days');
%     ylim([0 2]);
% 
%     subplot(2,3,5); hold on; % FA lick latency 
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     test_mr_fa = mr_fa(ctl,:);
%     test_mo_fa = mo_fa(ctl,:);
%     data = [test_mr_fa(:);test_mo_fa(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = test_mr_fa(i,:);
%         g2 = test_mo_fa(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     ylabel('FA lick latency (s)');xlabel('Days');
%     ylim([0 2]);
% 
%     subplot(2,3,6); hold on; % d primes
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mr_hit = rdprime(ctl,:);
%     ctl_mo_hit = odprime(ctl,:);
%     data = [ctl_mr_hit(:);ctl_mo_hit(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = ctl_mr_hit(i,:);
%         g2 = ctl_mo_hit(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylabel("d'");xlabel('Days');
%     ylim([-1 4]);
end

%% Figure for Kishore's R01


%%%%%%%%%%
% Let's correct d' in probe: fix after max
% pdrimeOrignial = pdprime;
% pdprime = pdrimeOrignial;
% [~,wm] = max(pdprime(:,1:10),[],2);
% for i=1:nSubj
%     pdprime(i,wm(i):end) = pdprime(i,wm(i));
% end    
%%%%%%%%%%

figure; hold on;
cd058=4;cd060=6;
ctl = find(explist==2);
% ctl(1) = []; % remove cd055
test = find(explist==1);
maxDays = 6;

% plot H and FA rate in PROBE of CD058 and CD060
subplot(2,3,1);hold on;
cd058_hit = subjrates{cd058}(1:maxDays,3)';
if ~isempty(find(diff(cd058_hit)<0))
    cd058_hit(find(diff(cd058_hit)<0)+1) = cd058_hit(find(diff(cd058_hit)<0,1));
end
cd058_fa = subjrates{cd058}(1:maxDays,4)';
cd060_hit = subjrates{cd060}(1:maxDays,3)';
% if ~isempty(find(diff(cd060_hit)<0))
%     cd060_hit(find(diff(cd060_hit)<0)+1) = cd060_hit(find(diff(cd060_hit)<0,1));
% end
cd060_fa = subjrates{cd060}(1:maxDays,4)';
plot(1:maxDays,cd058_hit,'ko');
% plot(1:nFiles,cd058_fa,'k--');
plot(1:maxDays,cd060_hit,'bo');

f = @(param,xval) param(1) + ( param(2)-param(1) )./ (   1 + 10.^( ( param(3) - xval ) * param(4) )   );
param = sigm_fit((1:maxDays)',cd058_hit',[],[],0);
x = (1:maxDays)';
y = cd058_hit';
x_vector=min(x):(max(x)-min(x))/100:max(x);
plot(x,y,'ko',x_vector,f(param,x_vector),'k-');

param = sigm_fit((1:maxDays)',cd060_hit',[],[],0);
y = cd060_hit';
plot(x,y,'bo',x_vector,f(param,x_vector),'b-');

% plot(1:nFiles,cd060_fa,'b--');
ylabel('Action proportion');xlabel('Days');
xlim([0.5 maxDays]);
   
% % plot d' in PROBE of CD058 and CD060
% DAY = 4;
% subplot(2,3,2);hold on;
% cd058_pdprime = pdprime(cd058,1:maxDays);
% if ~isempty(find(diff(cd058_pdprime)<0))
%     tochange = find(diff(cd058_pdprime)<0)+1;
%     tochange(tochange<3) = [];
%     cd058_pdprime(tochange) = cd058_pdprime(min(tochange)-1);
% end
% plot(1:maxDays,cd058_pdprime,'ko');
% plot(1:maxDays,pdprime(cd060,1:maxDays),'bo');
% 
% param = sigm_fit((1:maxDays)',cd058_pdprime',[],[],0);
% y = cd058_pdprime';
% plot(x,y,'ko',x_vector,f(param,x_vector),'k-');
% 
% param = sigm_fit((1:maxDays)',pdprime(cd060,1:maxDays)',[],[],0);
% y = pdprime(cd060,1:maxDays)';
% plot(x,y,'bo',x_vector,f(param,x_vector),'b-');
% 
% ylabel("d' probe");xlabel('Days');
% PlotIntervals([DAY-0.5 DAY+0.5]);
% xlim([0.5 maxDays]);
% ylim([-0.5 4]);
figure;
% plot d' in PROBE of TEST vs CTL
DAY = 5;
DAYS = 1:nFiles;
subplot(2,3,2);hold on;
% shadedErrorBar(DAYS,nanmean(pdprime(ctl,DAYS)),nansem(pdprime(ctl,DAYS)),{'color','k'},0);
% shadedErrorBar(DAYS,nanmean(pdprime(test,DAYS)),nansem(pdprime(test,DAYS)),{'color','b'},0);
f = @(param,xval) param(1) + ( param(2)-param(1) )./ (   1 + 10.^( ( param(3) - xval ) * param(4) )   );
x = DAYS';
for i=1:nSubj
    y = pdprime(i,DAYS);
    param = sigm_fit(x,y,[],[],0);
    x_vector = min(x):(max(x)-min(x))/100:max(x);
    if ismember(i,ctl)
        plot(x,y,'k.',x_vector,f(param,x_vector),'k-');
    else
        plot(x,y,'b.',x_vector,f(param,x_vector),'b-');
    end
end

% plot d' in PROBE of TEST vs CTL
subplot(2,3,3);hold on;
pdprime_ctl = pdprime(ctl,DAY);
pdprime_test = pdprime(test,DAY);
data = [pdprime_ctl(:);pdprime_test(:)];
groups = [ones(length(pdprime_ctl(:)),1);2*ones(length(pdprime_test(:)),1)];
barwitherrn(data,groups);hold on;
g1 = linspace(0.75,1.25,length(pdprime_ctl(:)));
g2 = linspace(1.75,2.25,length(pdprime_test(:)));
plot(g1,pdprime_ctl(:),'ko');
plot(g2,pdprime_test(:),'bo');
ylabel(["d' probe DAY" num2str(DAY)]);
xlim([0.5 2.5]);
set(gca,'xtick', [1 2],'xticklabel',{'CTL' 'TEST'});
title('nCtl=3; nTest=3');
% [~,p] = ttest2(data(groups==1),data(groups==2));
% [~,p] = ttest(data(groups==1),0);
% [~,p] = ttest(data(groups==2),0);

% % plot mean d' in REINFORCED ligh off + on for TEST and CTL
% subplot(2,2,3);hold on;
% shadedErrorBar(1:nFiles,nanmean(rfdprime(ctl,:)),nansem(rfdprime(ctl,:)),{'color','k'});
% shadedErrorBar(1:nFiles,nanmean(rfdprime(test,:)),nansem(rfdprime(test,:)),{'color','b'});
% xlim([1 maxDays]);
% ylabel("d' reinforced TEST vs CTL");xlabel('Days');
% title('nCtl=2; nTest=3');

% plot mean d' in REINFORCED Light ON for TEST vs CTL
DAYS = 1:nFiles;
subplot(2,4,5);hold on;
shadedErrorBar(DAYS,nanmean(odprime(ctl,DAYS)),nansem(odprime(ctl,DAYS)),{'color','k'},0);
shadedErrorBar(DAYS,nanmean(odprime(test,DAYS)),nansem(odprime(test,DAYS)),{'color','b'},0);
% plot(DAYS,nanmean(odprime(ctl,DAYS)),'k','linewidth',2);
% plot(DAYS,odprime(ctl,DAYS),'k');
% plot(DAYS,nanmean(odprime(test,DAYS)),'b','linewidth',2);
% plot(DAYS,odprime(test,DAYS),'b');
xlim([DAYS(1)-0.5 DAYS(end)]);
ylabel("d' reinforced Light ON, TEST vs CTL");xlabel('Days');
title('nCtl=3; nTest=3');

% plot mean d' in REINFORCED Light ON for TEST vs CTL
DAYS = 5:10;
% DAYS = 1:maxDays;
% colors = repelem(linspace(0,0.8,length(DAYS))',1,3);
subplot(2,4,6);hold on;
odprime_test = odprime(test,DAYS);
odprime_ctl = odprime(ctl,DAYS);
data = [odprime_ctl(:);odprime_test(:)];
groups = [ones(length(odprime_ctl(:)),1);2*ones(length(odprime_test(:)),1)];
barwitherrn(data,groups);hold on;
g1 = linspace(0.75,1.25,length(odprime_ctl(:)));
g2 = linspace(1.75,2.25,length(odprime_test(:)));
% for i=1:length(DAYS)
%     plot(g1(i), odprime_ctl(:,i), 'o-','color',colors(i,:));
%     plot(g2(i), odprime_test(:,i), 'o-','color',colors(i,:));
% end
plot(g1, odprime_ctl(:), 'ko');
plot(g2, odprime_test(:), 'ko');
[~,p] = ttest2(data(groups==1),data(groups==2));
% [~,p] = ttest(data(groups==1),0);
% [~,p] = ttest(data(groups==2),0);
set(gca,'xtick', [1 2],'xticklabel',{'CTL' 'TEST'});
ylabel("d' reinforced Light ON, TEST vs CTL");
xlabel(['Days ' num2str(DAYS(1)) ' to ' num2str(DAYS(end))]);
title(['nCtl=2; nTest=3, p=' num2str(p)]);

% plot mean d' in REINFORCED ligh off and in REINFORCED light on for TEST
DAYS = 1:nFiles;
subplot(2,4,7);hold on;
shadedErrorBar(DAYS,nanmean(rdprime(test,DAYS)),nansem(rdprime(test,DAYS)),{'color','k'},0);
shadedErrorBar(DAYS,nanmean(odprime(test,DAYS)),nansem(odprime(test,DAYS)),{'color','b'},0);
% plot(DAYS,nanmean(rdprime(test,DAYS)),'k','linewidth',2);
% plot(DAYS,rdprime(test,DAYS),'k');
% plot(DAYS,nanmean(odprime(test,DAYS)),'b','linewidth',2);
% plot(DAYS,odprime(test,DAYS),'b');
xlim([DAYS(1)-0.5 DAYS(end)]);
ylabel("d' reinforced for TEST");xlabel('Days');
title('Light OFF vs ON');

% plot mean d' in REINFORCED ligh off and in REINFORCED light on for TEST
%DAYS = 5:nFiles; % not sure why it starts from 5?
DAYS = 1:nFiles;
subplot(2,4,8);hold on;
rdprime_lightoff = rdprime(test,DAYS);
rdprime_lighton = odprime(test,DAYS);
data = [rdprime_lightoff(:);rdprime_lighton(:)];
groups = [ones(length(rdprime_lightoff(:)),1);2*ones(length(rdprime_lighton(:)),1)];
barwitherrn(data,groups);hold on;
% for i=1:length(DAYS)
%     plot(1:2, [rdprime_lightoff(:,i) rdprime_lighton(:,i)], 'o-','color',colors(i,:));
% end
plot(1:2, [rdprime_lightoff(:) rdprime_lighton(:)], 'ko-');
[~,p] = ttest(data(groups==1),data(groups==2));
set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
ylabel("d' reinforced for TEST");
xlabel(['Days ' num2str(DAYS(1)) ' to ' num2str(DAYS(end))]);
title(['Light OFF vs ON, p=' num2str(p)]);
%% Figure for lab meeting 3/6/2020

%%%%%%%%%%
% Let's correct d' in probe: fix after max
pdrimeOrignial = pdprime;
pdprime = pdrimeOrignial;
[~,wm] = max(pdprime(:,1:10),[],2);
for i=1:nSubj
    pdprime(i,wm(i):end) = pdprime(i,wm(i));
end    
%%%%%%%%%%

fig=figure; hold on;
ctl = find(explist==2);
test = find(explist==1);
maxDays = nFiles;

% plot H rate in PROBE
subplot(2,5,1);hold on;
for i=1:nSubj
    if ismember(i,ctl)
        plot(1:maxDays,subjrates{i}(1:maxDays,3),'k-');
    else
        plot(1:maxDays,subjrates{i}(1:maxDays,3),'b-');
    end
end
ylabel('HIT rate');xlabel('Days');
xlim([DAYS(1)-0.5 DAYS(end)]);
% plot FA rate in PROBE
subplot(2,5,2);hold on;
for i=1:nSubj
    if ismember(i,ctl)
        plot(1:maxDays,subjrates{i}(1:maxDays,4),'k-');
    else
        plot(1:maxDays,subjrates{i}(1:maxDays,4),'b-');
    end
end
ylabel('FA rate');xlabel('Days');
xlim([DAYS(1)-0.5 DAYS(end)]);
% plot dprime in PROBE
subplot(2,5,3);hold on;
plot(DAYS,pdrimeOrignial(ctl,:),'k');
plot(DAYS,pdrimeOrignial(test,:),'b');
ylim([-2 5]);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylabel("d'");xlabel('Days');
% plot dprime CORRECTED in PROBE
DAYS = 1:nFiles;
subplot(2,5,4);hold on;
% shadedErrorBar(DAYS,nanmean(pdprime(ctl,DAYS)),nansem(pdprime(ctl,DAYS)),{'color','k'},0);
% shadedErrorBar(DAYS,nanmean(pdprime(test,DAYS)),nansem(pdprime(test,DAYS)),{'color','b'},0);
f = @(param,xval) param(1) + ( param(2)-param(1) )./ (   1 + 10.^( ( param(3) - xval ) * param(4) )   );
x = DAYS';
for i=1:nSubj
    y = pdprime(i,DAYS);
    param = sigm_fit(x,y,[],[],0);
    x_vector = min(x):(max(x)-min(x))/100:max(x);
    if ismember(i,ctl)
        plot(x,y,'k.',x_vector,f(param,x_vector),'k-');
    else
        plot(x,y,'b.',x_vector,f(param,x_vector),'b-');
    end
end
ylim([-1 5]);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylabel("d' fixed");xlabel('Days');
DAY = 4;
PlotIntervals([DAY-0.5 DAY+0.5]);
% plot d' in PROBE of TEST vs CTL
subplot(2,5,5);hold on;
pdprime_ctl = pdprime(ctl,DAY);
pdprime_test = pdprime(test,DAY);
data = [pdprime_ctl(:);pdprime_test(:)];
groups = [ones(length(pdprime_ctl(:)),1);2*ones(length(pdprime_test(:)),1)];
barwitherrn(data,groups);hold on;
g1 = linspace(0.75,1.25,length(pdprime_ctl(:)));
g2 = linspace(1.75,2.25,length(pdprime_test(:)));
plot(g1,pdprime_ctl(:),'ko');
plot(g2,pdprime_test(:),'bo');
ylabel(["d' probe DAY" num2str(DAY)]);
xlim([0.5 2.5]);
set(gca,'xtick', [1 2],'xticklabel',{'CTL' 'TEST'});
title('nCtl=3; nTest=3');

% plot H rate in REINFORCED FULL
subplot(2,5,6);hold on;
for i=1:nSubj
    if ismember(i,ctl)
        plot(1:maxDays,subjrates{i}(1:maxDays,7),'k-');
    else
        plot(1:maxDays,subjrates{i}(1:maxDays,7),'b-');
    end
end
ylabel('HIT rate');xlabel('Days');
xlim([DAYS(1)-0.5 DAYS(end)]);
% plot FA rate in REINFORCED FULL
subplot(2,5,7);hold on;
for i=1:nSubj
    if ismember(i,ctl)
        plot(1:maxDays,subjrates{i}(1:maxDays,8),'k-');
    else
        plot(1:maxDays,subjrates{i}(1:maxDays,8),'b-');
    end
end
ylabel('FA rate');xlabel('Days');
xlim([DAYS(1)-0.5 DAYS(end)]);
% plot dprime in REINFORCED FULL
subplot(2,5,8);hold on;
plot(DAYS,rfdprime(ctl,:),'k');
plot(DAYS,rfdprime(test,:),'b');
ylim([-2 5]);
ylabel("d'");xlabel('Days');
% plot average dprime in REINFORCED FULL
subplot(2,5,9);hold on;
shadedErrorBar(DAYS,nanmean(rfdprime(ctl,:)),nansem(rfdprime(ctl,:)),{'color','k'},0.5);
shadedErrorBar(DAYS,nanmean(rfdprime(test,:)),nansem(rfdprime(test,:)),{'color','b'},0.5);
ylim([-2 5]);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylabel("d'");xlabel('Days');

savefig = false;
if savefig
    cd('T:\su\Figures\opto1\');
    saveas(fig,'probeVSreinforced.pdf');
end
%%

fig=figure; hold on; % REINFORCED LIGHT-ON vs LIGHT-OFF

subplot(2,4,1);hold on; % 1st row: CTL



DAY = 4;
PlotIntervals([DAY-0.5 DAY+0.5]);
% plot d' in PROBE of TEST vs CTL
subplot(2,5,5);hold on;
pdprime_ctl = pdprime(ctl,DAY);
pdprime_test = pdprime(test,DAY);
data = [pdprime_ctl(:);pdprime_test(:)];
groups = [ones(length(pdprime_ctl(:)),1);2*ones(length(pdprime_test(:)),1)];
barwitherrn(data,groups);hold on;
g1 = linspace(0.75,1.25,length(pdprime_ctl(:)));
g2 = linspace(1.75,2.25,length(pdprime_test(:)));
plot(g1,pdprime_ctl(:),'ko');
plot(g2,pdprime_test(:),'bo');
ylabel(["d' probe DAY" num2str(DAY)]);
xlim([0.5 2.5]);
set(gca,'xtick', [1 2],'xticklabel',{'CTL' 'TEST'});
title('nCtl=3; nTest=3');

[~,p] = ttest2(data(groups==1),data(groups==2));
p = ranksum(data(groups==1),data(groups==2));

% % plot mean d' in REINFORCED ligh off + on for TEST and CTL
% subplot(2,2,3);hold on;
% shadedErrorBar(1:nFiles,nanmean(rfdprime(ctl,:)),nansem(rfdprime(ctl,:)),{'color','k'});
% shadedErrorBar(1:nFiles,nanmean(rfdprime(test,:)),nansem(rfdprime(test,:)),{'color','b'});
% xlim([1 maxDays]);
% ylabel("d' reinforced TEST vs CTL");xlabel('Days');
% title('nCtl=2; nTest=3');

% plot mean d' in REINFORCED Light ON for TEST vs CTL
DAYS = 1:nFiles;
subplot(2,4,5);hold on;
shadedErrorBar(DAYS,nanmean(odprime(ctl,DAYS)),nansem(odprime(ctl,DAYS)),{'color','k'},0);
shadedErrorBar(DAYS,nanmean(odprime(test,DAYS)),nansem(odprime(test,DAYS)),{'color','b'},0);
% plot(DAYS,nanmean(odprime(ctl,DAYS)),'k','linewidth',2);
% plot(DAYS,odprime(ctl,DAYS),'k');
% plot(DAYS,nanmean(odprime(test,DAYS)),'b','linewidth',2);
% plot(DAYS,odprime(test,DAYS),'b');
xlim([DAYS(1)-0.5 DAYS(end)]);
ylabel("d' reinforced Light ON, TEST vs CTL");xlabel('Days');
title('nCtl=2; nTest=3');

% plot mean d' in REINFORCED Light ON for TEST vs CTL
DAYS = 5:10;
% DAYS = 1:maxDays;
% colors = repelem(linspace(0,0.8,length(DAYS))',1,3);
subplot(2,4,6);hold on;
odprime_test = odprime(test,DAYS);
odprime_ctl = odprime(ctl,DAYS);
data = [odprime_ctl(:);odprime_test(:)];
groups = [ones(length(odprime_ctl(:)),1);2*ones(length(odprime_test(:)),1)];
barwitherrn(data,groups);hold on;
g1 = linspace(0.75,1.25,length(odprime_ctl(:)));
g2 = linspace(1.75,2.25,length(odprime_test(:)));
% for i=1:length(DAYS)
%     plot(g1(i), odprime_ctl(:,i), 'o-','color',colors(i,:));
%     plot(g2(i), odprime_test(:,i), 'o-','color',colors(i,:));
% end
plot(g1, odprime_ctl(:), 'ko');
plot(g2, odprime_test(:), 'ko');
[~,p] = ttest2(data(groups==1),data(groups==2));
% [~,p] = ttest(data(groups==1),0);
% [~,p] = ttest(data(groups==2),0);
set(gca,'xtick', [1 2],'xticklabel',{'CTL' 'TEST'});
ylabel("d' reinforced Light ON, TEST vs CTL");
xlabel(['Days ' num2str(DAYS(1)) ' to ' num2str(DAYS(end))]);
title(['nCtl=2; nTest=3, p=' num2str(p)]);

% plot mean d' in REINFORCED ligh off and in REINFORCED light on for TEST
DAYS = 1:nFiles;
subplot(2,4,7);hold on;
shadedErrorBar(DAYS,nanmean(rdprime(test,DAYS)),nansem(rdprime(test,DAYS)),{'color','k'},0);
shadedErrorBar(DAYS,nanmean(odprime(test,DAYS)),nansem(odprime(test,DAYS)),{'color','b'},0);
% plot(DAYS,nanmean(rdprime(test,DAYS)),'k','linewidth',2);
% plot(DAYS,rdprime(test,DAYS),'k');
% plot(DAYS,nanmean(odprime(test,DAYS)),'b','linewidth',2);
% plot(DAYS,odprime(test,DAYS),'b');
xlim([DAYS(1)-0.5 DAYS(end)]);
ylabel("d' reinforced for TEST");xlabel('Days');
title('Light OFF vs ON');

% plot mean d' in REINFORCED ligh off and in REINFORCED light on for TEST
DAYS = 5:nFiles;
subplot(2,4,8);hold on;
rdprime_lightoff = rdprime(test,DAYS);
rdprime_lighton = odprime(test,DAYS);
data = [rdprime_lightoff(:);rdprime_lighton(:)];
groups = [ones(length(rdprime_lightoff(:)),1);2*ones(length(rdprime_lighton(:)),1)];
barwitherrn(data,groups);hold on;
% for i=1:length(DAYS)
%     plot(1:2, [rdprime_lightoff(:,i) rdprime_lighton(:,i)], 'o-','color',colors(i,:));
% end
plot(1:2, [rdprime_lightoff(:) rdprime_lighton(:)], 'ko-');
[~,p] = ttest(data(groups==1),data(groups==2));
set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
ylabel("d' reinforced for TEST");
xlabel(['Days ' num2str(DAYS(1)) ' to ' num2str(DAYS(end))]);
title(['Light OFF vs ON, p=' num2str(p)]);


%% Averages: action rate, d', % correct

% %%%%%%%%%
% Let's correct d' in probe: fix after max
% pdrimeOrignial = pdprime;
% pdprime = pdrimeOrignial;
% [~,wm] = max(pdprime(:,1:13),[],2);
% for i=1:nSubj
%     pdprime(i,wm(i):end) = pdprime(i,wm(i));
% end    
% %%%%%%%%%

% explist = [2 1 2 2 1 1 2 1 2 1 1 2]; % 1 = test, 2 = ctl
% cond = 'all';
% explist = [2 1 2 2 1 1 2 1 0 1 0 2]; % only no learner removed
% explist = [2 1 2 2 1 1 2 1 0 1 0 2 1 1 1 2]; % 1 = test, 2 = ctl
% cond = 'nonlearnerremoved';
% explist = [2 1 2 2 1 0 2 1 0 1 0 0 1 1 1 2]; % all bad removed
% cond = 'badremoved';
% cd060 (test): good then bad
% cd065 (test): no learning at all
% cd063 (ctl): no learning, only one day with perf>1
% cd066 (ctl): good then bad

% explist = [2 1 2 1 1 2 1 2 2 1 2 1 1 1 2 1]; % 90-10, cohorts 1+2
% explist = [0 0 0 1 0 2 1 2 2 1 2 1 1 1 2 1]; % 90-10, cohorts 1+2, non learner removed
% explist = [0 0 0 1 0 2 1 0 2 1 2 1 1 1 2 1]; % 90-10, cohorts 1+2, non learner + non stable removed

% subjlist = {'CD055','CD057','CD058','CD061','CD063','CD066','CD076','CD077','CD079','CD082','CD084','CD085','CD087','CD091'};
% explist = [1 1 1 1 0 1 1 0 0 2 2 2 2 2]; % 50-50 + 90-10, control only, non learner removed

% cond = 'nonlearnerremoved';
cond = 'badremoved';
ctl = find(explist==2);
test = find(explist==1);
savefig = false;

fig=figure; hold on;
DAYS = 1:17;

% plot HIT and FA in REINFORCED ligh off and in REINFORCED light on for
% TEST and CTL
subplot(4,4,1);hold on;title('Test');
test_m_r = nanmean(test_subjrates_fa_r);
test_s_r = nansem(test_subjrates_fa_r);
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'--','color','k'},0.5);
test_m_o = nanmean(test_subjrates_fa_o);
test_s_o = nansem(test_subjrates_fa_o);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'--','color','b'},0.5);

% figure; 
% plot(DAYS,test_subjrates_fa_r(:,DAYS),'k--');hold on;
% plot(DAYS,test_subjrates_hit_r(:,DAYS),'k');

test_m_r = nanmean(test_subjrates_hit_r);
test_s_r = nansem(test_subjrates_hit_r);
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','k'},0.5);
test_m_o = nanmean(test_subjrates_hit_o);
test_s_o = nansem(test_subjrates_hit_o);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'-','color','b'},0.5);
ylim([0 1]);

subplot(4,4,2);hold on;title('Ctl');
ctl_m_r = nanmean(ctl_subjrates_fa_r);
ctl_s_r = nansem(ctl_subjrates_fa_r);
shadedErrorBar(DAYS,ctl_m_r(DAYS),ctl_s_r(DAYS),{'--','color','k'},0.5);
ctl_m_o = nanmean(ctl_subjrates_fa_o);
ctl_s_o = nansem(ctl_subjrates_fa_o);
shadedErrorBar(DAYS,ctl_m_o(DAYS),ctl_s_o(DAYS),{'--','color','b'},0.5);

ctl_m_r = nanmean(ctl_subjrates_hit_r);
ctl_s_r = nansem(ctl_subjrates_hit_r);
shadedErrorBar(DAYS,ctl_m_r(DAYS),ctl_s_r(DAYS),{'-','color','k'},0.5);
ctl_m_o = nanmean(ctl_subjrates_hit_o);
ctl_s_o = nansem(ctl_subjrates_hit_o);
shadedErrorBar(DAYS,ctl_m_o(DAYS),ctl_s_o(DAYS),{'-','color','b'},0.5);
ylim([0 1]); 

subplot(4,4,3);hold on;title('Light OFF');
test_m_r = nanmean(test_subjrates_fa_r);
test_s_r = nansem(test_subjrates_fa_r);
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'--','color','b'},0.5);
test_m_r = nanmean(test_subjrates_hit_r);
test_s_r = nansem(test_subjrates_hit_r);
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','b'},0.5);

ctl_m_r = nanmean(ctl_subjrates_fa_r);
ctl_s_r = nansem(ctl_subjrates_fa_r);
shadedErrorBar(DAYS,ctl_m_r(DAYS),ctl_s_r(DAYS),{'--','color','k'},0.5);
ctl_m_r = nanmean(ctl_subjrates_hit_r);
ctl_s_r = nansem(ctl_subjrates_hit_r);
shadedErrorBar(DAYS,ctl_m_r(DAYS),ctl_s_r(DAYS),{'-','color','k'},0.5);
ylim([0 1]);

subplot(4,4,4);hold on;title('Light ON');
test_m_o = nanmean(test_subjrates_fa_o);
test_s_o = nansem(test_subjrates_fa_o);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'--','color','b'},0.5);
test_m_o = nanmean(test_subjrates_hit_o);
test_s_o = nansem(test_subjrates_hit_o);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'-','color','b'},0.5);

ctl_m_o = nanmean(ctl_subjrates_fa_o);
ctl_s_o = nansem(ctl_subjrates_fa_o);
shadedErrorBar(DAYS,ctl_m_o(DAYS),ctl_s_o(DAYS),{'--','color','k'},0.5);
ctl_m_o = nanmean(ctl_subjrates_hit_o);
ctl_s_o = nansem(ctl_subjrates_hit_o);
shadedErrorBar(DAYS,ctl_m_o(DAYS),ctl_s_o(DAYS),{'-','color','k'},0.5);
ylim([0 1]);

subplot(4,4,5);hold on;
shadedErrorBar(DAYS,nanmean(rdprime(test,DAYS)),nansem(rdprime(test,DAYS)),{'color','k'},0.5);
shadedErrorBar(DAYS,nanmean(odprime(test,DAYS)),nansem(odprime(test,DAYS)),{'color','b'},0.5);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylim([0 4.5]);

subplot(4,4,6);hold on;
shadedErrorBar(DAYS,nanmean(rdprime(ctl,DAYS)),nansem(rdprime(ctl,DAYS)),{'color','k'},0.5);
shadedErrorBar(DAYS,nanmean(odprime(ctl,DAYS)),nansem(odprime(ctl,DAYS)),{'color','b'},0.5);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylim([0 4.5]);

subplot(4,4,7);hold on;
shadedErrorBar(DAYS,nanmean(rdprime(ctl,DAYS)),nansem(rdprime(ctl,DAYS)),{'color','k'},0.5);
shadedErrorBar(DAYS,nanmean(rdprime(test,DAYS)),nansem(rdprime(test,DAYS)),{'color','b'},0.5);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylim([0 4.5]);

subplot(4,4,8);hold on;
shadedErrorBar(DAYS,nanmean(odprime(ctl,DAYS)),nansem(odprime(ctl,DAYS)),{'color','k'},0.5);
shadedErrorBar(DAYS,nanmean(odprime(test,DAYS)),nansem(odprime(test,DAYS)),{'color','b'},0.5);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylim([0 4.5]);
title(['Light ON, ctl (n=' num2str(length(ctl)) ', black) VS test (n=' num2str(length(test)) ', blue)'])
subplot(4,4,9);hold on;title('Test');
test_m_r = nanmean(((1-test_subjrates_fa_r)+test_subjrates_hit_r)/2);
test_s_r = nansem(((1-test_subjrates_fa_r)+test_subjrates_hit_r));
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','k'},0.5);
test_m_o = nanmean(((1-test_subjrates_fa_o)+test_subjrates_hit_o)/2);
test_s_o = nansem(((1-test_subjrates_fa_o)+test_subjrates_hit_o)/2);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'-','color','b'},0.5);
ylim([0.4 1]);

subplot(4,4,10);hold on;title('Ctl');
test_m_r = nanmean(((1-ctl_subjrates_fa_r)+ctl_subjrates_hit_r)/2);
test_s_r = nansem(((1-ctl_subjrates_fa_r)+ctl_subjrates_hit_r));
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','k'},0.5);
test_m_o = nanmean(((1-ctl_subjrates_fa_o)+ctl_subjrates_hit_o)/2);
test_s_o = nansem(((1-ctl_subjrates_fa_o)+ctl_subjrates_hit_o)/2);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'-','color','b'},0.5);
ylim([0.4 1]);

subplot(4,4,11);hold on;title('Light OFF');
test_m_r = nanmean(((1-ctl_subjrates_fa_r)+ctl_subjrates_hit_r)/2);
test_s_r = nansem(((1-ctl_subjrates_fa_r)+ctl_subjrates_hit_r));
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','k'},0.5);
test_m_r = nanmean(((1-test_subjrates_fa_r)+test_subjrates_hit_r)/2);
test_s_r = nansem(((1-test_subjrates_fa_r)+test_subjrates_hit_r));
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','b'},0.5);
ylim([0.4 1]);

subplot(4,4,12);hold on;title('Light ON');
test_m_o = nanmean(((1-ctl_subjrates_fa_o)+ctl_subjrates_hit_o)/2);
test_s_o = nansem(((1-ctl_subjrates_fa_o)+ctl_subjrates_hit_o)/2);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'-','color','k'},0.5);
ylim([0.4 1]);
test_m_o = nanmean(((1-test_subjrates_fa_o)+test_subjrates_hit_o)/2);
test_s_o = nansem(((1-test_subjrates_fa_o)+test_subjrates_hit_o)/2);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'-','color','b'},0.5);
ylim([0.4 1]);


ctl_subjrates_hit_p = nan(length(ctl),nFiles);
ctl_subjrates_fa_p = nan(length(ctl),nFiles);
test_subjrates_hit_p = nan(length(test),nFiles);
test_subjrates_fa_p = nan(length(test),nFiles);
c = 1 ; t = 1;
for i=1:nSubj    
    if ismember(i,ctl)
        nFiles = length(subjrates{i}(:,3)');
        ctl_subjrates_hit_p(c,1:nFiles) = subjrates{i}(:,3)';
        ctl_subjrates_fa_p(c,1:nFiles) = subjrates{i}(:,4)';
        c=c+1;
    elseif ismember(i,test)
        nFiles = length(subjrates{i}(:,3)');
        test_subjrates_hit_p(t,1:nFiles) = subjrates{i}(:,3)';
        test_subjrates_fa_p(t,1:nFiles) = subjrates{i}(:,4)';
        t=t+1;
    end       
end
nFiles = maxFiles;
subplot(4,3,10);hold on;title('Test vs Ctl, probe');
test_m_r = nanmean(test_subjrates_fa_p);
test_s_r = nansem(test_subjrates_fa_p);
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'--','color','b'},0.5);
test_m_o = nanmean(ctl_subjrates_fa_p);
test_s_o = nansem(ctl_subjrates_fa_p);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'--','color','k'},0.5);

test_m_r = nanmean(test_subjrates_hit_p);
test_s_r = nansem(test_subjrates_hit_p);
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','b'},0.5);
test_m_o = nanmean(ctl_subjrates_hit_p);
test_s_o = nansem(ctl_subjrates_hit_p);
shadedErrorBar(DAYS,test_m_o(DAYS),test_s_o(DAYS),{'-','color','k'},0.5);
ylim([0 1]);

subplot(4,3,11);hold on; title('dprime probe');
shadedErrorBar(DAYS,nanmean(pdprime(ctl,DAYS)),nansem(pdprime(ctl,DAYS)),{'color','k'},0.5);
shadedErrorBar(DAYS,nanmean(pdprime(test,DAYS)),nansem(pdprime(test,DAYS)),{'color','b'},0.5);
xlim([DAYS(1)-0.5 DAYS(end)]);
ylim([0 4.5]);

subplot(4,3,12);hold on; title('% correct probe');
test_m_r = nanmean(((1-ctl_subjrates_fa_p)+ctl_subjrates_hit_p)/2);
test_s_r = nansem(((1-ctl_subjrates_fa_p)+ctl_subjrates_hit_p));
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','k'},0.5);
test_m_r = nanmean(((1-test_subjrates_fa_p)+test_subjrates_hit_p)/2);
test_s_r = nansem(((1-test_subjrates_fa_p)+test_subjrates_hit_p));
shadedErrorBar(DAYS,test_m_r(DAYS),test_s_r(DAYS),{'-','color','b'},0.5);
ylim([0.4 1]);

if savefig
    cd(pathsave);
    saveas(fig,['Average-Action-Performance_Correct-Reinforced-Probe-' cond '.pdf']);
    close(fig);
end

%% KISHORE'S GRANT MCNIGHT
% Averages: action rate, d', % correct (export svg compatible)

% %%%%%%%%%
% Let's correct d' in probe: fix after max
% pdrimeOrignial = pdprime;
% pdprime = pdrimeOrignial;
% [~,wm] = max(pdprime(:,1:13),[],2);
% for i=1:nSubj
%     pdprime(i,wm(i):end) = pdprime(i,wm(i));
% end    
% %%%%%%%%%

% explistfull = [2 1 2 1 1 2 1 2 2 1 2 1 1 1 2 1]; % 90-10, cohorts 1+2
% explist = [0 0 0 1 0 2 1 2 2 1 2 1 1 1 2 1]; % 90-10, cohorts 1+2, non learner removed
explist = [0 0 0 1 0 2 1 0 2 1 2 1 1 1 2 1]; % 90-10, cohorts 1+2, non learner + non stable removed

% subjlist = {'CD055','CD057','CD058','CD061','CD063','CD066','CD076','CD077','CD079','CD082','CD084','CD085','CD087','CD091'};
% explist = [1 1 1 1 0 1 1 0 0 2 2 2 2 2]; % 50-50 + 90-10, control only, non learner removed

% cond = 'nonlearnerremoved';
cond = 'badremoved';
ctl = find(explist==2);
test = find(explist==1);
savefig = false;

fig=figure; hold on;
DAYS = 1:21;

% plot HIT and FA in REINFORCED ligh off and in REINFORCED light on for
% TEST and CTL
subplot(4,4,1);hold on;title('Test, FA');
test_m_r = nanmean(test_subjrates_fa_r);
test_s_r = nansem(test_subjrates_fa_r);
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),test_m_r(DAYS)+test_s_r(DAYS),'color','k');
test_m_o = nanmean(test_subjrates_fa_o);
test_s_o = nansem(test_subjrates_fa_o);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),test_m_o(DAYS)+test_s_o(DAYS),'color','b');
ylim([0 1]);

subplot(4,4,2);hold on;title('Test, H');
test_m_r = nanmean(test_subjrates_hit_r);
test_s_r = nansem(test_subjrates_hit_r);
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),test_m_r(DAYS)+test_s_r(DAYS),'color','k');
test_m_o = nanmean(test_subjrates_hit_o);
test_s_o = nansem(test_subjrates_hit_o);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),test_m_o(DAYS)+test_s_o(DAYS),'color','b');
ylim([0 1]);

subplot(4,4,3);hold on;title('Ctl, FA');
ctl_m_r = nanmean(ctl_subjrates_fa_r);
ctl_s_r = nansem(ctl_subjrates_fa_r);
PlotMean(DAYS,ctl_m_r(DAYS),ctl_m_r(DAYS)-ctl_s_r(DAYS),ctl_m_r(DAYS)+ctl_s_r(DAYS),'color','k');
ctl_m_o = nanmean(ctl_subjrates_fa_o);
ctl_s_o = nansem(ctl_subjrates_fa_o);
PlotMean(DAYS,ctl_m_o(DAYS),ctl_m_o(DAYS)-ctl_s_o(DAYS),ctl_m_o(DAYS)+ctl_s_o(DAYS),'color','b');
ylim([0 1]); 

subplot(4,4,4);hold on;title('Ctl, H');
ctl_m_r = nanmean(ctl_subjrates_hit_r);
ctl_s_r = nansem(ctl_subjrates_hit_r);
PlotMean(DAYS,ctl_m_r(DAYS),ctl_m_r(DAYS)-ctl_s_r(DAYS),ctl_m_r(DAYS)+ctl_s_r(DAYS),'color','k');
ctl_m_o = nanmean(ctl_subjrates_hit_o);
ctl_s_o = nansem(ctl_subjrates_hit_o);
PlotMean(DAYS,ctl_m_o(DAYS),ctl_m_o(DAYS)-ctl_s_o(DAYS),ctl_m_o(DAYS)+ctl_s_o(DAYS),'color','b');
ylim([0 1]); 

subplot(4,4,5);hold on;title('Light OFF, Test');
test_m_r = nanmean(test_subjrates_fa_r);
test_s_r = nansem(test_subjrates_fa_r);
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),test_m_r(DAYS)+test_s_r(DAYS),'color','b');
test_m_r = nanmean(test_subjrates_hit_r);
test_s_r = nansem(test_subjrates_hit_r);
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),test_m_r(DAYS)+test_s_r(DAYS),'color','b');

subplot(4,4,6);hold on;title('Light OFF, Ctl');
ctl_m_r = nanmean(ctl_subjrates_fa_r);
ctl_s_r = nansem(ctl_subjrates_fa_r);
PlotMean(DAYS,ctl_m_r(DAYS),ctl_m_r(DAYS)-ctl_s_r(DAYS),ctl_m_r(DAYS)+ctl_s_r(DAYS),'color','k');
ctl_m_r = nanmean(ctl_subjrates_hit_r);
ctl_s_r = nansem(ctl_subjrates_hit_r);
PlotMean(DAYS,ctl_m_r(DAYS),ctl_m_r(DAYS)-ctl_s_r(DAYS),ctl_m_r(DAYS)+ctl_s_r(DAYS),'color','k');
ylim([0 1]);

subplot(4,4,7);hold on;title('Light ON, Test');
test_m_o = nanmean(test_subjrates_fa_o);
test_s_o = nansem(test_subjrates_fa_o);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),test_m_o(DAYS)+test_s_o(DAYS),'color','b');
test_m_o = nanmean(test_subjrates_hit_o);
test_s_o = nansem(test_subjrates_hit_o);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),test_m_o(DAYS)+test_s_o(DAYS),'color','b');
ylim([0 1]);

subplot(4,4,8);hold on;title('Light ON, Ctl');
ctl_m_o = nanmean(ctl_subjrates_fa_o);
ctl_s_o = nansem(ctl_subjrates_fa_o);
PlotMean(DAYS,ctl_m_o(DAYS),ctl_m_o(DAYS)-ctl_s_o(DAYS),ctl_m_o(DAYS)+ctl_s_o(DAYS),'color','k');
ctl_m_o = nanmean(ctl_subjrates_hit_o);
ctl_s_o = nansem(ctl_subjrates_hit_o);
PlotMean(DAYS,ctl_m_o(DAYS),ctl_m_o(DAYS)-ctl_s_o(DAYS),ctl_m_o(DAYS)+ctl_s_o(DAYS),'color','k');
ylim([0 1]);

subplot(4,4,9);hold on;title("Light OFF, d'");
PlotMean(DAYS,nanmean(rdprime(test,DAYS)),...
    nanmean(rdprime(test,DAYS))-nansem(rdprime(test,DAYS)),...
    nanmean(rdprime(test,DAYS))+nansem(rdprime(test,DAYS)),...
    'color','b');
PlotMean(DAYS,nanmean(rdprime(ctl,DAYS)),...
    nanmean(rdprime(ctl,DAYS))-nansem(rdprime(ctl,DAYS)),...
    nanmean(rdprime(ctl,DAYS))+nansem(rdprime(ctl,DAYS)),...
    'color','k');
ylim([0 4.5]);
xlim([0 DAYS(end)]);

subplot(4,4,10);hold on;title("Light ON, d'");
PlotMean(DAYS,nanmean(odprime(test,DAYS)),...
    nanmean(odprime(test,DAYS))-nansem(odprime(test,DAYS)),...
    nanmean(odprime(test,DAYS))+nansem(odprime(test,DAYS)),...
    'color','b');
PlotMean(DAYS,nanmean(odprime(ctl,DAYS)),...
    nanmean(odprime(ctl,DAYS))-nansem(odprime(ctl,DAYS)),...
    nanmean(odprime(ctl,DAYS))+nansem(odprime(ctl,DAYS)),...
    'color','k');
ylim([0 4.5]);
xlim([0 DAYS(end)]);

subplot(4,4,11);hold on;title('Light OFF, %');
test_m_r = nanmean(((1-ctl_subjrates_fa_r)+ctl_subjrates_hit_r)/2);
test_s_r = nansem(((1-ctl_subjrates_fa_r)+ctl_subjrates_hit_r));
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),...
    test_m_r(DAYS)+test_s_r(DAYS),'color','k');
test_m_r = nanmean(((1-test_subjrates_fa_r)+test_subjrates_hit_r)/2);
test_s_r = nansem(((1-test_subjrates_fa_r)+test_subjrates_hit_r));
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),...
    test_m_r(DAYS)+test_s_r(DAYS),'color','b');
ylim([0.4 1]);
xlim([0 DAYS(end)]);

subplot(4,4,12);hold on;title('Light ON, %');
test_m_o = nanmean(((1-ctl_subjrates_fa_o)+ctl_subjrates_hit_o)/2);
test_s_o = nansem(((1-ctl_subjrates_fa_o)+ctl_subjrates_hit_o)/2);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),...
    test_m_o(DAYS)+test_s_o(DAYS),'color','k');
ylim([0.4 1]);
test_m_o = nanmean(((1-test_subjrates_fa_o)+test_subjrates_hit_o)/2);
test_s_o = nansem(((1-test_subjrates_fa_o)+test_subjrates_hit_o)/2);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),...
    test_m_o(DAYS)+test_s_o(DAYS),'color','b');
ylim([0.4 1]);
xlim([0 DAYS(end)]);

ctl_subjrates_hit_p = nan(length(ctl),nFiles);
ctl_subjrates_fa_p = nan(length(ctl),nFiles);
test_subjrates_hit_p = nan(length(test),nFiles);
test_subjrates_fa_p = nan(length(test),nFiles);
c = 1 ; t = 1;
for i=1:nSubj    
    if ismember(i,ctl)
        nFiles = length(subjrates{i}(:,3)');
        ctl_subjrates_hit_p(c,1:nFiles) = subjrates{i}(:,3)';
        ctl_subjrates_fa_p(c,1:nFiles) = subjrates{i}(:,4)';
        c=c+1;
    elseif ismember(i,test)
        nFiles = length(subjrates{i}(:,3)');
        test_subjrates_hit_p(t,1:nFiles) = subjrates{i}(:,3)';
        test_subjrates_fa_p(t,1:nFiles) = subjrates{i}(:,4)';
        t=t+1;
    end       
end
maxdayprobe =6;
nFiles = maxFiles;
subplot(4,4,13);hold on;title('Test-Ctl, FA probe');
test_m_r = nanmean(test_subjrates_fa_p);
test_s_r = nansem(test_subjrates_fa_p);
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),...
    test_m_r(DAYS)+test_s_r(DAYS),'color','b');
test_m_o = nanmean(ctl_subjrates_fa_p);
test_s_o = nansem(ctl_subjrates_fa_p);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),...
    test_m_o(DAYS)+test_s_o(DAYS),'color','k');
xlim([1 maxdayprobe]);
ylim([0 1]);

subplot(4,4,14);hold on;title('Test-Ctl, H probe');
test_m_r = nanmean(test_subjrates_hit_p);
test_s_r = nansem(test_subjrates_hit_p);
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),...
    test_m_r(DAYS)+test_s_r(DAYS),'color','b');
test_m_o = nanmean(ctl_subjrates_hit_p);
test_s_o = nansem(ctl_subjrates_hit_p);
PlotMean(DAYS,test_m_o(DAYS),test_m_o(DAYS)-test_s_o(DAYS),...
    test_m_o(DAYS)+test_s_o(DAYS),'color','k');
ylim([0 1]);
xlim([1 maxdayprobe]);

subplot(4,4,15);hold on; title('dprime probe');
PlotMean(DAYS,nanmean(pdprime(ctl,DAYS)),...
    nanmean(pdprime(ctl,DAYS))-nansem(pdprime(ctl,DAYS)),...
    nanmean(pdprime(ctl,DAYS))+nansem(pdprime(ctl,DAYS)),'color','k');
PlotMean(DAYS,nanmean(pdprime(test,DAYS)),...
    nanmean(pdprime(test,DAYS))-nansem(pdprime(test,DAYS)),...
    nanmean(pdprime(test,DAYS))+nansem(pdprime(test,DAYS)),'color','b');
xlim([1 maxdayprobe]);
ylim([0 3.5]);

subplot(4,4,16);hold on; title('% correct probe');
test_m_r = nanmean(((1-ctl_subjrates_fa_p)+ctl_subjrates_hit_p)/2);
test_s_r = nansem(((1-ctl_subjrates_fa_p)+ctl_subjrates_hit_p));
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),...
    test_m_r(DAYS)+test_s_r(DAYS),'color','k');
test_m_r = nanmean(((1-test_subjrates_fa_p)+test_subjrates_hit_p)/2);
test_s_r = nansem(((1-test_subjrates_fa_p)+test_subjrates_hit_p));
PlotMean(DAYS,test_m_r(DAYS),test_m_r(DAYS)-test_s_r(DAYS),...
    test_m_r(DAYS)+test_s_r(DAYS),'color','b');
ylim([0.4 1]);
xlim([1 maxdayprobe]);


% if savefig
%     cd(pathsave);
%     saveas(fig,['Average-Action-Performance_Correct-Reinforced-Probe-' cond '.pdf']);
%     close(fig);
% end

%% Plot indiv action rate in reinforced
for i=1:size(ctl_subjrates_fa_r,1)
    figure;hold on;
    plot(DAYS,ctl_subjrates_fa_r(i,DAYS),'k--');
    plot(DAYS,ctl_subjrates_fa_o(i,DAYS),'b--');
    plot(DAYS,ctl_subjrates_hit_r(i,DAYS),'k');
    plot(DAYS,ctl_subjrates_hit_o(i,DAYS),'b');
    ylim([0 1]);
end
%% Plot indiv action rate in reinforced
c=0;t=0;
f = @(param,xval) param(1) + ( param(2)-param(1) )./ (   1 + 10.^( ( param(3) - xval ) * param(4) )   );
DAYS = 1:6;
x = DAYS;
colorhit = [47 196 182]/255;
colorfa = [255 159 25]/255;
for i=1:nSubj
    if ismember(i,ctl)
        figure;hold on;
        c = c+1;
        plot(DAYS,ctl_subjrates_fa_p(c,DAYS),'color',colorfa);  
%         plot(DAYS,ctl_subjrates_fa_p(c,DAYS),'o','color',colorfa);  
%         param = sigm_fit(x,ctl_subjrates_fa_p(c,DAYS),[],[],0);
%         fittedsig = f(param,x);
%         plot(x,fittedsig,'-','color',colorfa);
              
        plot(DAYS,ctl_subjrates_hit_p(c,DAYS),'color',colorhit);
%         plot(DAYS,ctl_subjrates_hit_p(c,DAYS),'o','color',colorhit); 
%         param = sigm_fit(x,ctl_subjrates_hit_p(c,DAYS),[],[],0);
%         fittedsig = f(param,x);
%         plot(x,fittedsig,'-','color',colorhit);
        
        ylim([0 1]);xlim([0 max(x)+1]);
        title([subjlist{i} ' - ' expnames{2}]);
    elseif ismember(i,test)
        figure;hold on;
        t = t+1;
        plot(DAYS,test_subjrates_fa_p(t,DAYS),'color',colorfa);
%         plot(DAYS,test_subjrates_fa_p(t,DAYS),'o','color',colorfa);        
%         param = sigm_fit(x,test_subjrates_fa_p(t,DAYS),[],[],0);
%         fittedsig = f(param,x);
%         plot(x,fittedsig,'-','color',colorfa);
                
        plot(DAYS,test_subjrates_hit_p(t,DAYS),'color',colorhit);    
%         plot(DAYS,test_subjrates_hit_p(t,DAYS),'o','color',colorhit);        
%         param = sigm_fit(x,test_subjrates_hit_p(t,DAYS),[],[],0);
%         fittedsig = f(param,x);
%         plot(x,fittedsig,'-','color',colorhit);
        
        ylim([0 1]);xlim([0 max(x)+1]);
        title([subjlist{i} ' - ' expnames{1}]);
    end
end

%% Look at difference in acquisition closer (probe)

savefig = false;

fig=figure; hold on;
days = 1:21;

% day of maximum d' on raw data 
[~,wm] = max(pdrimeOrignial,[],2);
daymaxperf = wm([ctl';test']);
groups = [ones(length(ctl),1);2*ones(length(test),1)];
subplot(2,4,1);hold on;
anovabar(daymaxperf,groups);hold on;
plot(groups,daymaxperf,'k.');
ylabel("Max raw d' (day)");
set(gca,'xtick',[1 2],'xticklabel',{'CTL','TEST'});

% day of maximum % correct on raw data
pc_ctl = ((1-ctl_subjrates_fa_p)+ctl_subjrates_hit_p)/2;
pc_test = ((1-test_subjrates_fa_p)+test_subjrates_hit_p)/2;
[~,daymaxperf] = max([pc_ctl;pc_test],[],2);
groups = [ones(length(ctl),1);2*ones(length(test),1)];
subplot(2,4,2);hold on;
anovabar(daymaxperf,groups);hold on;
plot(groups,daymaxperf,'k.');
ylabel('Max raw %correct (day)');
set(gca,'xtick',[1 2],'xticklabel',{'CTL','TEST'});

% day of maximum corrected d'
[~,wm] = max(pdprime,[],2);
daymaxperf = wm([ctl';test']);
groups = [ones(length(ctl),1);2*ones(length(test),1)];
subplot(2,4,3);hold on;
anovabar(daymaxperf,groups);hold on;
plot(groups,daymaxperf,'k.');
ylim([0 20]);
ylabel("Max corrected d' (day)");
set(gca,'xtick',[1 2],'xticklabel',{'CTL','TEST'});

% day of maximum perf on sigmoid fit on 
f = @(param,xval) param(1) + ( param(2)-param(1) )./ (   1 + 10.^( ( param(3) - xval ) * param(4) )   );
x = days';
y = pdprime([ctl';test'],days)';
x_vector = x;
fittedsig = nan(size(y,2),length(x_vector));
% subplot(2,4,4);hold on;cla
for i=1:size(y,2)
    param = sigm_fit(x,y(:,i),[],[],0);
%     if groups(i)==1
%         plot(x,y(:,i),'k.',x_vector,f(param,x_vector),'k-');
%     else
%         plot(x,y(:,i),'b.',x_vector,f(param,x_vector),'b-');
%     end
    fittedsig(i,:) = f(param,x_vector);
end
[~,daymaxperf] = max(round(fittedsig,1),[],2);
groups = [ones(length(ctl),1);2*ones(length(test),1)];
subplot(2,4,4);hold on;
anovabar(daymaxperf,groups);hold on;
plot(groups,daymaxperf,'k.');
ylim([0 20]);
ylabel("Max fitted d' (day)");
set(gca,'xtick',[1 2],'xticklabel',{'CTL','TEST'});

% day where raw d' >1
wm = nan(nSubj,1);
for i=1:nSubj
    if sum(pdrimeOrignial(i,:)'>1)~=0
        wm(i) = find(pdrimeOrignial(i,:)'>1,1);
    end
end
daymaxperf = wm([ctl';test']);
groups = [ones(length(ctl),1);2*ones(length(test),1)];
subplot(2,4,5);hold on;
anovabar(daymaxperf,groups);hold on;
plot(groups,daymaxperf,'k.');
ylim([0 20]);

wm = nan(nSubj,1);
for i=1:nSubj
    if sum(pdrimeOrignial(i,:)'>1.5)~=0
        wm(i) = find(pdrimeOrignial(i,:)'>1.5,1);
    end
end
daymaxperf = wm([ctl';test']);
groups = [ones(length(ctl),1);2*ones(length(test),1)];
anovabar(daymaxperf,groups+2);hold on;
plot(groups+2,daymaxperf,'k.');

wm = nan(nSubj,1);
for i=1:nSubj
    fd = find(pdrimeOrignial(i,:)'>2,1);
    if ~isempty(fd)
        wm(i) = find(pdrimeOrignial(i,:)'>2,1);
    end
end
daymaxperf = wm([ctl';test']);
groups = [ones(length(ctl),1);2*ones(length(test),1)];
anovabar(daymaxperf,groups+4);hold on;
plot(groups+4,daymaxperf,'k.');
ylim([0 10]);
ylabel("Day where d' > 1,1.5,2");

% corrected d' diff according to days
subplot(2,4,6);hold on;
for d=1:7
    daymaxperf = pdprime([ctl';test'],d);
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
end
ylim([0 3.5]);
ylabel("Corrected d', ctl and test");
xlabel('Days');

% %correct according to days
pc = [pc_ctl;pc_test];
subplot(2,4,7);hold on;
for d=1:7
    daymaxperf = pc(:,d);
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
    [~,p] = ttest2(daymaxperf(groups==1),daymaxperf(groups==2));
    disp(['day ' num2str(d) ', p=' num2str(p)]);
end
ylim([0.3 1]);
ylabel('% Correct, ctl and test');
xlabel('Days');

% corrected d' bloc1 diff according to days
% %%%%%%%%%%%%
% % Let's correct d' in bloc 1 of probe: fix after max
% pdrimeOrignial1 = pdprime1;
% pdprime1 = pdrimeOrignial1;
% [~,wm] = max(pdprime1(:,1:10),[],2);
% for i=1:nSubj
%     pdprime1(i,wm(i):end) = pdprime1(i,wm(i));
% end    
% % Let's correct d' in bloc 2 of probe: fix after max
% pdrimeOrignial2 = pdprime2;
% pdprime2 = pdrimeOrignial2;
% [~,wm] = max(pdprime2(:,1:10),[],2);
% for i=1:nSubj
%     pdprime2(i,wm(i):end) = pdprime2(i,wm(i));
% end    
% %%%%%%%%%%%%%
prevgroups = groups;

subplot(2,4,8);hold on;
for d=1:7
    daymaxperf = [pdprime1([ctl';test'],d) pdprime2([ctl';test'],d)]';
    daymaxperf = daymaxperf(:);    
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(1:2:end-1),'g.');hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(2:2:end),'k.');hold on;
end
ylim([0 2.5]);
ylabel("d', ctl and test, bloc 1&2");
xlabel('Days');

if savefig
    cd(pathsave);
    saveas(fig,['Average-ProbeDiff-CtlTest-1-' cond '.pdf']);
    close(fig);
end

fig=figure;hold on;
subplot(2,2,1);hold on;cla
for d=1:7
    daymaxperf = pdprime1([ctl';test'],d);
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'g.');hold on;
end
ylim([0 2.5]);
ylabel("d', ctl and test, bloc 1");
xlabel('Days');

subplot(2,2,2);hold on;cla
for d=1:7
    daymaxperf = pdprime2([ctl';test'],d);
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
end
ylim([0 2.5]);
ylabel("d', ctl and test, bloc 2");
xlabel('Days');

subplot(2,2,3);hold on;cla
for d=1:7
    daymaxperf = [pdprime1([ctl';test'],d) pdprime2([ctl';test'],d)]';
    daymaxperf = daymaxperf(:);    
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(1:2:end-1),'g.');hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(2:2:end),'k.');hold on;
end
ylim([0 2.5]);
ylabel("d', ctl and test, bloc 1&2");
xlabel('Days');

subplot(2,2,4);hold on;cla
for d=1:7
    daymaxperf = pdprime([ctl';test'],d);
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
end
ylim([0 3.5]);
ylabel("corrected d', ctl and test");
xlabel('Days');

if savefig
    cd(pathsave);
    saveas(fig,['Average-ProbeDiff-CtlTest-DprimeBloc1&2-' cond '.pdf']);
    close(fig);
end

ctl_subjrates_phit1 = nan(length(ctl),nFiles);
ctl_subjrates_phit2 = nan(length(ctl),nFiles);
ctl_subjrates_pfa1 = nan(length(ctl),nFiles);
ctl_subjrates_pfa2 = nan(length(ctl),nFiles);

test_subjrates_phit1 = nan(length(test),nFiles);
test_subjrates_phit2 = nan(length(test),nFiles);
test_subjrates_pfa1 = nan(length(test),nFiles);
test_subjrates_pfa2 = nan(length(test),nFiles);

c = 1 ; t = 1;
for i=1:nSubj    
    if ismember(i,ctl)
        nFiles = length(subjrates{i}(:,1)');
        ctl_subjrates_phit1(c,1:nFiles) = subjrates{i}(:,11)';
        ctl_subjrates_phit2(c,1:nFiles) = subjrates{i}(:,12)';
        
        ctl_subjrates_pfa1(c,1:nFiles) = subjrates{i}(:,13)';
        ctl_subjrates_pfa2(c,1:nFiles) = subjrates{i}(:,14)';
        c=c+1;
    elseif ismember(i,test)
        nFiles = length(subjrates{i}(:,1)');
        test_subjrates_phit1(t,1:nFiles) = subjrates{i}(:,11)';
        test_subjrates_phit2(t,1:nFiles) = subjrates{i}(:,12)';
        
        test_subjrates_pfa1(t,1:nFiles) = subjrates{i}(:,13)';
        test_subjrates_pfa2(t,1:nFiles) = subjrates{i}(:,14)';
        t=t+1;
    end       
end
nFiles = maxFiles;

pc_ctl_1 = ((1-ctl_subjrates_pfa1)+ctl_subjrates_phit1)/2;
pc_ctl_2 = ((1-ctl_subjrates_pfa2)+ctl_subjrates_phit2)/2;

pc_test_1 = ((1-test_subjrates_pfa1)+test_subjrates_phit1)/2;
pc_test_2 = ((1-test_subjrates_pfa2)+test_subjrates_phit2)/2;

fig=figure;hold on;
subplot(2,2,1);hold on;cla
for d=1:7
    daymaxperf = [pc_ctl_1(:,d);pc_test_1(:,d)];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'g.');hold on;
end
ylim([0.2 1]);
ylabel('% correct, ctl and test, bloc 1');
xlabel('Days');

subplot(2,2,2);hold on;cla
for d=1:7
    daymaxperf = [pc_ctl_2(:,d);pc_test_2(:,d)];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
end
ylim([0.2 1]);
ylabel('% correct, ctl and test, bloc 2');
xlabel('Days');

subplot(2,2,3);hold on;cla
for d=1:7
    daymaxperf = [[pc_ctl_1(:,d);pc_test_1(:,d)]  [pc_ctl_2(:,d);pc_test_2(:,d)]]';
    daymaxperf = daymaxperf(:);    
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(1:2:end-1),'g.');hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(2:2:end),'k.');hold on;
end
ylim([0.2 1]);
ylabel('% correct, ctl and test, bloc 1&2');
xlabel('Days');

subplot(2,2,4);hold on;cla
pc = [pc_ctl;pc_test];
for d=1:7
    daymaxperf = pc(:,d);
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
    [~,p] = ttest2(daymaxperf(groups==1),daymaxperf(groups==2));
    disp(['day ' num2str(d) ', p=' num2str(p)]);
end
ylim([0.2 1]);
ylabel('% correct, ctl and test, combined');
xlabel('Days');

if savefig
    cd(pathsave);
    saveas(fig,['Average-ProbeDiff-CtlTest-PercentCorrectBloc1&2-' cond '.pdf']);
    close(fig);
end

% ACTION RATE DIFFERENCES 
pc_ctl_1 = ctl_subjrates_phit1 - ctl_subjrates_pfa1;
pc_ctl_2 = ctl_subjrates_phit2 - ctl_subjrates_pfa2;

pc_test_1 = test_subjrates_phit1 - test_subjrates_pfa1;
pc_test_2 = test_subjrates_phit2 - test_subjrates_pfa2;

fig = figure;hold on;
subplot(2,2,1);hold on;
for d=1:7
    daymaxperf = [pc_ctl_1(:,d);pc_test_1(:,d)];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'g.');hold on;
end
ylim([-1 1]);
ylabel('Action rate difference (hit-fa), ctl and test, bloc 1');
xlabel('Days');
subplot(2,2,2);hold on;
for d=1:7
    daymaxperf = [pc_ctl_2(:,d);pc_test_2(:,d)];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
end
ylim([-1 1]);
ylabel('Action rate difference (hit-fa), ctl and test, bloc 2');
xlabel('Days');
subplot(2,2,3);hold on;
for d=1:7
    daymaxperf = [[pc_ctl_1(:,d);pc_test_1(:,d)]  [pc_ctl_2(:,d);pc_test_2(:,d)]]';
    daymaxperf = daymaxperf(:);    
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(1:2:end-1),'g.');hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(2:2:end),'k.');hold on;
end
ylim([-1 1]);
ylabel('Action rate difference (hit-fa), ctl and test, bloc 1&2');
xlabel('Days');

if savefig
    cd(pathsave);
    saveas(fig,['Average-ProbeDiff-CtlTest-DifferenceActionRateBloc1&2-' cond '.pdf']);
    close(fig);
end

% Days grouped two-by-two

% PERCENTAGES
pc_ctl_1 = ((1-ctl_subjrates_pfa1)+ctl_subjrates_phit1)/2;
pc_ctl_2 = ((1-ctl_subjrates_pfa2)+ctl_subjrates_phit2)/2;

pc_test_1 = ((1-test_subjrates_pfa1)+test_subjrates_phit1)/2;
pc_test_2 = ((1-test_subjrates_pfa2)+test_subjrates_phit2)/2;
ylimm = [0.2 1];

fig=figure;hold on;
subplot(2,2,1);hold on;cla
for d=1:2:6
    daymaxperfctl = pc_ctl_1(:,[d d+1]);
    daymaxperftest = pc_test_1(:,[d d+1]);
    daymaxperf = [daymaxperfctl(:);daymaxperftest(:)];
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'g.');hold on;
end
ylim(ylimm);
ylabel('%correct, ctl and test, bloc 1');
xlabel('Days grouped by two');

subplot(2,2,2);hold on;cla
for d=1:2:6
    daymaxperfctl = pc_ctl_2(:,[d d+1]);
    daymaxperftest = pc_test_2(:,[d d+1]);
    daymaxperf = [daymaxperfctl(:);daymaxperftest(:)];
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
end
ylim(ylimm);
ylabel('%correct, ctl and test, bloc 2');
xlabel('Days grouped by two');

prevgroups = groups;
subplot(2,2,3);hold on;cla
for d=1:2:6
    daymaxperfctl1 = pc_ctl_1(:,[d d+1]);
    daymaxperfctl2 = pc_ctl_2(:,[d d+1]);
    daymaxperfctl = [daymaxperfctl1(:);daymaxperfctl2(:)];
    daymaxperftest1 = pc_test_1(:,[d d+1]);
    daymaxperftest2 = pc_test_2(:,[d d+1]);
    daymaxperftest = [daymaxperftest1(:);daymaxperftest2(:)];    
    daymaxperf = [daymaxperfctl;daymaxperftest];
    groups = [ones(length(ctl)*2*2,1);2*ones(length(test)*2*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(1:2:end-1),'g.');hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(2:2:end),'k.');hold on;
    [~,p] = ttest2(daymaxperf(groups==1),daymaxperf(groups==2));
    disp(['day ' num2str(d) ' and ' num2str(d+1) ', p = ' num2str(p)]);
end
ylim(ylimm);
ylabel('%correct, ctl and test, bloc 1&2');
xlabel('Days grouped by two');

if savefig
    cd(pathsave);
    saveas(fig,['Average-ProbeDiff-CtlTest-PercentCorrectBloc1&2-DaysGrouped-' cond '.pdf']);
    close(fig);
end

% ACTION RATE DIFFERENCES
pc_ctl_1 = ctl_subjrates_phit1 - ctl_subjrates_pfa1;
pc_ctl_2 = ctl_subjrates_phit2 - ctl_subjrates_pfa2;

pc_test_1 = test_subjrates_phit1 - test_subjrates_pfa1;
pc_test_2 = test_subjrates_phit2 - test_subjrates_pfa2;
ylimm = [-0.5 1];

fig=figure;hold on;
subplot(2,2,1);hold on;cla
for d=1:2:6
    daymaxperfctl = pc_ctl_1(:,[d d+1]);
    daymaxperftest = pc_test_1(:,[d d+1]);
    daymaxperf = [daymaxperfctl(:);daymaxperftest(:)];
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'g.');hold on;
end
ylim(ylimm);
ylabel('Action rate difference (hit-fa), ctl and test, bloc 1');
xlabel('Days grouped by two');

subplot(2,2,2);hold on;cla
for d=1:2:6
    daymaxperfctl = pc_ctl_2(:,[d d+1]);
    daymaxperftest = pc_test_2(:,[d d+1]);
    daymaxperf = [daymaxperfctl(:);daymaxperftest(:)];
    groups = [ones(length(ctl)*2,1);2*ones(length(test)*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,daymaxperf,'k.');hold on;
end
ylim(ylimm);
ylabel('Action rate difference (hit-fa), ctl and test, bloc 2');
xlabel('Days grouped by two');

prevgroups = groups;
subplot(2,2,3);hold on;cla
for d=1:2:6
    daymaxperfctl1 = pc_ctl_1(:,[d d+1]);
    daymaxperfctl2 = pc_ctl_2(:,[d d+1]);
    daymaxperfctl = [daymaxperfctl1(:);daymaxperfctl2(:)];
    daymaxperftest1 = pc_test_1(:,[d d+1]);
    daymaxperftest2 = pc_test_2(:,[d d+1]);
    daymaxperftest = [daymaxperftest1(:);daymaxperftest2(:)];    
    daymaxperf = [daymaxperfctl;daymaxperftest];
    groups = [ones(length(ctl)*2*2,1);2*ones(length(test)*2*2,1)];
    anovabar(daymaxperf,groups+(d-1)*2);hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(1:2:end-1),'g.');hold on;
    plot(prevgroups+(d-1)*2,daymaxperf(2:2:end),'k.');hold on;
    [~,p] = ttest2(daymaxperf(groups==1),daymaxperf(groups==2));
    disp(['day ' num2str(d) ' and ' num2str(d+1) ', p = ' num2str(p)]);
end
ylim(ylimm);
title(['ctl n=' num2str(sum(groups==1)) ', test n=' num2str(sum(groups==2))]);
ylabel('Action rate difference (hit-fa), ctl and test, bloc 1&2');
xlabel('Days grouped by two');

if savefig
    cd(pathsave);
    saveas(fig,['Average-ProbeDiff-CtlTest-DifferenceActionRateBloc1&2-DaysGrouped-' cond '.pdf']);
    close(fig);
end

%% Look at lick rate and latency in PROBE
% explist = [2 1 2 2 1 1 2 1 2 1 1 2]; % 1 = test, 2 = ctl
% cond = 'all';
% explist = [2 1 2 2 1 1 2 1 0 1 0 2]; % only no learner removed
% cond = 'nonlearnerremoved';
% explist = [2 1 2 2 1 0 2 1 0 1 0 0]; % all bad removed
% cond = 'badremoved';
ctl = find(explist==2);
test = find(explist==1);
savefig = false;

maxday = 17;
days = 1:maxday;
ndays = length(days);

fig=figure;
subplot(2,4,1);hold on;
shadedErrorBar(1:nFiles,nanmean(mp_hit(ctl,:)),nansem(mp_hit(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mp_hit(test,:)),nansem(mp_hit(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 2.5]);
xlabel('Days');
ylabel('HIT Lick latency (s)');
subplot(2,4,2);hold on;
for d=1:ndays
    datactl = mp_hit(ctl,d);
    datatest = mp_hit(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

subplot(2,4,3);hold on;
shadedErrorBar(1:nFiles,nanmean(mp_fa(ctl,:)),nansem(mp_fa(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mp_fa(test,:)),nansem(mp_fa(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 2.5]);
xlabel('Days');
ylabel('FA Lick latency (s)');
subplot(2,4,4);hold on;
for d=1:ndays
    datactl = mp_fa(ctl,d);
    datatest = mp_fa(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

subplot(2,4,5);hold on;
shadedErrorBar(1:nFiles,nanmean(mlp_hit(ctl,:)),nansem(mlp_hit(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mlp_hit(test,:)),nansem(mlp_hit(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 4]);
xlabel('Days');
ylabel('HIT Lick rate (Hz)');
subplot(2,4,6);hold on;
for d=1:ndays
    datactl = mlp_hit(ctl,d);
    datatest = mlp_hit(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

subplot(2,4,7);hold on;
shadedErrorBar(1:nFiles,nanmean(mlp_fa(ctl,:)),nansem(mlp_fa(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mlp_fa(test,:)),nansem(mlp_fa(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 4]);
xlabel('Days');
ylabel('FA Lick rate (Hz)');
subplot(2,4,8);hold on;
for d=1:ndays
    datactl = mlp_fa(ctl,d);
    datatest = mlp_fa(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

if savefig
    cd(pathsave);
    saveas(fig,['Average-Probe-LickLatency-LickRate-' cond '.pdf']);
    close(fig);
end

%% Look at lick rate and latency in REINFORCED
% explist = [2 1 2 2 1 1 2 1 2 1 1 2]; % 1 = test, 2 = ctl
% cond = 'all';
% explist = [2 1 2 2 1 1 2 1 0 1 0 2]; % only no learner removed
% cond = 'nonlearnerremoved';
% explist = [2 1 2 2 1 0 2 1 0 1 0 0]; % all bad removed
% cond = 'badremoved';
ctl = find(explist==2);
test = find(explist==1);
savefig = false;

maxday = 21;
days = 1:maxday;
ndays = length(days);

fig=figure;
subplot(2,4,1);hold on;
% shadedErrorBar(1:nFiles,nanmean(mr_hit(ctl,:)),nansem(mr_hit(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mr_hit(test,:)),nansem(mr_hit(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 2.5]);
xlabel('Days');
ylabel('HIT Lick latency (s)');
subplot(2,4,2);hold on;
for d=1:ndays
    datactl = mr_hit(ctl,d);
    datatest = mr_hit(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

subplot(2,4,3);hold on;
shadedErrorBar(1:nFiles,nanmean(mr_fa(ctl,:)),nansem(mr_fa(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mr_fa(test,:)),nansem(mr_fa(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 2.5]);
xlabel('Days');
ylabel('FA Lick latency (s)');
subplot(2,4,4);hold on;
for d=1:ndays
    datactl = mr_fa(ctl,d);
    datatest = mr_fa(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

subplot(2,4,5);hold on;
shadedErrorBar(1:nFiles,nanmean(mlr_hit(ctl,:)),nansem(mlr_hit(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mlr_hit(test,:)),nansem(mlr_hit(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 7]);
xlabel('Days');
ylabel('HIT Lick rate (Hz)');
subplot(2,4,6);hold on;
for d=1:ndays
    datactl = mlr_hit(ctl,d);
    datatest = mlr_hit(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

subplot(2,4,7);hold on;
shadedErrorBar(1:nFiles,nanmean(mlr_fa(ctl,:)),nansem(mlr_fa(ctl,:)),{'k'},0.3);
shadedErrorBar(1:nFiles,nanmean(mlr_fa(test,:)),nansem(mlr_fa(test,:)),{'b'},0.3);
xlim([1 maxday]);ylim([0 7]);
xlabel('Days');
ylabel('FA Lick rate (Hz)');
subplot(2,4,8);hold on;
for d=1:ndays
    datactl = mlr_fa(ctl,d);
    datatest = mlr_fa(test,d);
    data = [datactl;datatest];
    groups = [ones(length(ctl),1);2*ones(length(test),1)];
    anovabar(data,groups+(d-1)*2);hold on;
    plot(groups+(d-1)*2,data,'k.');hold on;
end

if savefig
    cd(pathsave);
    saveas(fig,['Average-ReinforcedLightOFF-LickLatency-LickRate-' cond '.pdf']);
    close(fig);
end


%% Averages lick PSTHs

% explist = [2 1 2 2 1 1 2 1 2 1 1 2]; % 1 = test, 2 = ctl
% cond = 'all';
% explist = [2 1 2 2 1 1 2 1 0 1 0 2]; % only no learner removed
% cond = 'nonlearnerremoved';
% explist = [2 1 2 2 1 0 2 1 0 1 0 0]; % all bad removed
% cond = 'badremoved';
ctl = find(explist==2);
test = find(explist==1);
savefig = false;

smoothfactor = 0;
maxday = 21;
l_r_hit_ctl=[];
l_r_fa_ctl=[];
l_o_hit_ctl=[];
l_o_fa_ctl=[];
l_p_hit_ctl=[];
l_p_fa_ctl=[];
l_r_hit_test=[];
l_r_fa_test=[];
l_o_hit_test=[];
l_o_fa_test=[];
l_p_hit_test=[];
l_p_fa_test=[];
groups_ctl=[];
groups_test=[];
for i=1:nSubj
    if ismember(i,ctl)
        l_r_hit_ctl = [l_r_hit_ctl;lickhistrhit{i}];
        l_r_fa_ctl = [l_r_fa_ctl;lickhistrfa{i}];
        l_o_hit_ctl = [l_o_hit_ctl;lickhistohit{i}];
        l_o_fa_ctl = [l_o_fa_ctl;lickhistofa{i}];
        l_p_hit_ctl = [l_p_hit_ctl;lickhistphit{i}];
        l_p_fa_ctl = [l_p_fa_ctl;lickhistpfa{i}];
        groups_ctl = [groups_ctl;(1:size(lickhistpfa{i},1))'];
    elseif ismember(i,test)
        l_r_hit_test = [l_r_hit_test;lickhistrhit{i}];
        l_r_fa_test = [l_r_fa_test;lickhistrfa{i}];
        l_o_hit_test = [l_o_hit_test;lickhistohit{i}];
        l_o_fa_test = [l_o_fa_test;lickhistofa{i}];
        l_p_hit_test = [l_p_hit_test;lickhistphit{i}];
        l_p_fa_test = [l_p_fa_test;lickhistpfa{i}];
        groups_test = [groups_test;(1:size(lickhistpfa{i},1))'];
    end
end
x = linspace(-1,4,nbins);

ylimm_hit = [0 80];
ylimm_fa = [0 80];
ylimm_phit = [0 6];
ylimm_pfa = [0 6];

fig=figure; % reinfHIT, ctl and test
for d=1:maxday
    subplot(5,5,d);hold on;
    shadedErrorBar(x,Smooth(nanmean(l_r_hit_ctl(groups_ctl==d,:)),smoothfactor),Smooth(nansem(l_r_hit_ctl(groups_ctl==d,:)),smoothfactor),{'k'},0.3);
    shadedErrorBar(x,Smooth(nanmean(l_r_hit_test(groups_test==d,:)),smoothfactor),Smooth(nansem(l_r_hit_test(groups_test==d,:)),smoothfactor),{'b'},0.3);
    ylim(ylimm_hit);
end
if savefig
    cd(pathsave);
    saveas(fig,['Average-LickPSTH-reinforced-hit-' cond '.pdf']);
    close(fig);
end

fig=figure; % reinfFA, ctl and test
for d=1:maxday
    subplot(5,5,d);hold on;
    shadedErrorBar(x,Smooth(mean(l_r_fa_ctl(groups_ctl==d,:)),smoothfactor),Smooth(sem(l_r_fa_ctl(groups_ctl==d,:)),smoothfactor),{'k--'},0.3);
    shadedErrorBar(x,Smooth(mean(l_r_fa_test(groups_test==d,:)),smoothfactor),Smooth(sem(l_r_fa_test(groups_test==d,:)),smoothfactor),{'b--'},0.3);
    ylim(ylimm_fa);
end
if savefig
    cd(pathsave);
    saveas(fig,['Average-LickPSTH-reinforced-fa-' cond '.pdf']);
    close(fig);
end

x = linspace(-1,4,length(bins));
fig=figure; % optoHIT, ctl and test
for d=1:maxday
    subplot(5,5,d);hold on;
    shadedErrorBar(x,Smooth(nanmean(l_o_hit_ctl(groups_ctl==d,:)),smoothfactor),Smooth(nansem(l_o_hit_ctl(groups_ctl==d,:)),smoothfactor),{'k'},0.3);
    shadedErrorBar(x,Smooth(nanmean(l_o_hit_test(groups_test==d,:)),smoothfactor),Smooth(nansem(l_o_hit_test(groups_test==d,:)),smoothfactor),{'b'},0.3);
    ylim(ylimm_hit);
end
if savefig
    cd(pathsave);
    saveas(fig,['Average-LickPSTH-opto-hit-' cond '.pdf']);
    close(fig);
end

fig=figure; % optoFA, ctl and test
for d=1:maxday
    subplot(5,5,d);hold on;
    shadedErrorBar(x,Smooth(mean(l_o_fa_ctl(groups_ctl==d,:)),smoothfactor),Smooth(sem(l_o_fa_ctl(groups_ctl==d,:)),smoothfactor),{'k--'},0.3);
    shadedErrorBar(x,Smooth(mean(l_o_fa_test(groups_test==d,:)),smoothfactor),Smooth(sem(l_o_fa_test(groups_test==d,:)),smoothfactor),{'b--'},0.3);
    ylim(ylimm_fa);
end
if savefig
    cd(pathsave);
    saveas(fig,['Average-LickPSTH-opto-fa-' cond '.pdf']);
    close(fig);
end

fig=figure; % probeHIT, ctl and test
for d=1:maxday
    subplot(5,5,d);hold on;
    shadedErrorBar(x,Smooth(nanmean(l_p_hit_ctl(groups_ctl==d,:)),smoothfactor),Smooth(nansem(l_p_hit_ctl(groups_ctl==d,:)),smoothfactor),{'k'},0.3);
    shadedErrorBar(x,Smooth(nanmean(l_p_hit_test(groups_test==d,:)),smoothfactor),Smooth(nansem(l_p_hit_test(groups_test==d,:)),smoothfactor),{'b'},0.3);
    ylim(ylimm_phit);
end
if savefig
    cd(pathsave);
    saveas(fig,['Average-LickPSTH-probe-hit-' cond '.pdf']);
    close(fig);
end

fig=figure; % probeFA, ctl and test
for d=1:maxday
    subplot(5,5,d);hold on;
    shadedErrorBar(x,Smooth(mean(l_p_fa_ctl(groups_ctl==d,:)),smoothfactor),Smooth(sem(l_p_fa_ctl(groups_ctl==d,:)),smoothfactor),{'k--'},0.3);
    shadedErrorBar(x,Smooth(mean(l_p_fa_test(groups_test==d,:)),smoothfactor),Smooth(sem(l_p_fa_test(groups_test==d,:)),smoothfactor),{'b--'},0.3);
    ylim(ylimm_pfa);
end
if savefig
    cd(pathsave);
    saveas(fig,['Average-LickPSTH-probe-fa-' cond '.pdf']);
    close(fig);
end

%% Look at trial history

SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
H=1;M=2;FA=3;CR=4;
DATA = cell(nSubj,2);

%%
for nbSubj = 1:nSubj
    
    matrix = MAT{nbSubj};
    % e.g. Probability to have a HIT after a HIT

    % days = [1;5;10;15;20];
    days = 1:max(matrix(:,SESS));
    nDays = length(days);
    outcs = [H;M;FA;CR]; 
    outcs_tone = [1 2;1 2;3 4;3 4];
    data = nan(length(outcs),length(outcs),nDays);
    outcNames = {'HIT';'MISS';'FA';'CR'};
    for d=1:nDays

        sess1 = matrix(matrix(:,SESS)==days(d) & matrix(:,CTXT)==2,:);
        nTrials = size(sess1,1);
        nPrevs = 1;
        idx = repmat((1:nPrevs+1)',1,nTrials) + repmat((0:nTrials-1),nPrevs+1,1);
        idx = idx';
        idx(idx(:,end)>nTrials,:) = [];
        idx_outcome = sess1(idx,OUTCOME);
        idx_outcome = reshape(idx_outcome,size(idx,1),size(idx,2));

    %     figure;
        for i=1:length(outcs)
    %         subplot(2,2,i);hold on;
            outc = outcs(i);
            outc_tone = outcs_tone(i,:);
            outc_given_prev = nan(length(outcs),1);
            for j=1:length(outcs)
                prev = outcs(j);
                outc_given_prev(j,:) = sum( idx_outcome(idx_outcome(:,1:end-1)==prev,end)==outc ) / sum(idx_outcome(idx_outcome(:,1:end-1)==prev,end)==outc_tone(1) | idx_outcome(idx_outcome(:,1:end-1)==prev,end)==outc_tone(2));
            end
    %         plot(outc_given_prev);
    %         ylabel(['Proba ' outcNames{i}]);
    %         set(gca,'xtick',[1;2;3;4],'xticklabels',{'H','M','FA','CR'});
    %         ylim([0 1]);
    %         if i==1, title([subjlist{nbSubj} ' - Day ' num2str(days(d))]); end
            data(i,:,d) = outc_given_prev';
        end
    end
    % figure;
    % colors = {'g','b','r','k'};
    % for i=1:length(outcs)
    %     subplot(2,2,i);
    %     toplot = squeeze(data(i,:,:));
    %     for j=1:4
    %         plot(toplot(j,:),'color',colors{j},'DisplayName',outcNames{j});hold on;
    %     end
    %     ylabel(['Proba ' outcNames{i}]);
    %     xlabel('Days');
    %     ylim([0 1]);
    %     if i==1, title(subjlist{nbSubj}); end
    % end
    disp(num2str(nbSubj));
    DATA{nbSubj,2} = data;

end
%% 
figure;
colors = {'g','b','r','k'};
for i=1:length(outcs)
    subplot(2,2,i);
    toplot = squeeze(data(i,:,:));
    for j=1:4
        plot(toplot(j,:),'color',colors{j},'DisplayName',outcNames{j});hold on;
    end
    ylabel(['Proba ' outcNames{i}]);
    xlabel('Days');
    ylim([0 1]);
    if i==1, title(subjlist{nbSubj}); end
end

%% Clean data
