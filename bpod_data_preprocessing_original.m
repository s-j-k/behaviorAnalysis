%To pre-process bpod GNG data, requires all sessions for each mouse in a mouse
%identifator folder and all mice folders within one general
%"folder'=folder. Additionally, put a string that contains all animal
%numbers to be analyzed within that folder. Outputs data as a structure that
%contains all relevant information per mouse with all mice sessions in
%cells. Creates a folder within the general folder, called "Analysis" where
%a .mat version of all data is saved.
%S.M. 11/25/2019 created
%S.M. 03/20/2020 updated needs to correct the fact that dA2 does not allow different number of files to be bound  
%S.M. 04/23/2020 corrected latest version but needs to correct that if within one day there were different sessions whit different tones played,
%it only saves the tone on the first session
%S.M. 09/03/2020 corrected and now it saves the tones that were played in
%each session

clear all
addpath('\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\Bpod_scripts') %path of your preprocessing functions
% addpath('/Users/shar/Documents/APP_project')
% addpath('/Users/shar/Documents/Analysis_APP_project/Behavior_bpod')
% addpath('A:\Sharlen\Behavior\CitricAcid')
addpath('A:\Sharlen\Behavior\CitricAcid') %path where your data is stored
% addpath('\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\APP_Project\Behavior')
% addpath('\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\Behavior\DSP4')
% addpath('\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\Behavior\astro_DREADDs')
% addpath('A:\Sharlen\Behavior\CitricAcid\bpodraw')

%AARON'S MICE: 
% AW={'AW017', 'AW018', 'AW019', 'AW020', 'AW021', 'AW022','AW023', 'AW024', 'AW025', 'AW026', 'AW027', 'AW028', 'AW029', 'AW030', 'AW031'}; 
% %KELLY'S MICE: 
% KF={'KF022', 'KF023', 'KF024', 'KF026', 'KF028'};
% %FANGCHEN's MICE: 
% FZ={'FZ007', 'FZ008', 'FZ009', 'FZ010'}; %'FZ001', 'FZ002', 'FZ003', 'FZ004', 'FZ005', 'FZ006', 
% %AM MICE: 
% AM={'AM015', 'AM017', 'AM018', 'AM019', 'AM020', 'AM022', 'AM023', 'AM024'};

%CITRIC ACID VS WATER RESTRICTION MICE COHORT 2 (after covid):
% CA={'SM039', 'SM040', 'SM041', 'SM042', 'SM043', 'SM044', 'SM045', 'SM046'};
% dob={'04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020'};
% genotype={'WR','WR','CA','CA','WR','WR','CA', 'CA'};

%CITRIC ACID VS WATER RESTRICTION MICE COHORT 3 (after covid):
% CA2={'SM049', 'SM050', 'SM051', 'SM052', 'SM053', 'SM054', 'SM055', 'SM056'};
% dob={'06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020'};
% genotype={'CA','CA','CA','CA','WR','WR','WR','WR'};

%APP MICE 10-12MO COHORT 
% APP10to12={'SM062', 'SM063', 'SM064', 'SM065', 'SM066', 'SM067'};
% dob={'01/07/2020','01/07/2020','01/07/2020','01/07/2020','01/07/2020','01/16/2020'};
% sex=[2,2,2,2,2,2];
% genotype=[2,2,1,1,1,1];    

%APP MICE 6-8MO COHORT 
% APP6to8={'SM068', 'SM069', 'SM070'};
% dob={'04/26/2020','04/26/2020','04/30/2020'};
% sex=[2,2,2];
% genotype=[1,2,1];    

%DSP4 vs Saline mice (Zyan/Angel)
% DSP4={'AL001', 'AL002', 'ZW001','ZW002','ZW003'};
% dob={'10/11/2020','10/11/2020','10/11/2020','10/11/2020','10/11/2020'};
% sex=[2,2,2,2,2]; %2 male, 1 female
% genotype=[1,2,1,2,2]; %1=saline injected, 2=dsp4 injected
% rawdata_folder='\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\Behavior\DSP4\';
% GNG_protocol='LickingGNG_Normal_cam'; %'LickingGNG_Reverse_cam';
% save_name='data_GNG_Normal'; %'data_GNG_Reverse';
% save_folder=rawdata_folder;

%CNO vs Saline mice (GFAP-dreadds Gi)
% AstroDREADDs={'SM072', 'SM073', 'SM075','SM076','SM077'};
% dob={'10/11/2020','10/11/2020','10/11/2020','10/11/2020','10/11/2020'};
% sex=[2,2,2,2,2]; %2 male, 1 female
% genotype=[1,2,2,1,2]; %1=saline injected, 2=CNO injected
% rawdata_folder='\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\Behavior\astro_DREADDs\';
% GNG_protocol='LickingGNG_Normal_cam';
% save_name='data_GNGNormal';
% save_folder=rawdata_folder;

CA2={'SM079', 'SM080', 'SM081','SM082','SM083', 'SM084', 'SM085', 'SM086','SM087','SM088','SM089','SM090'};
dob={'02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021','02/16/2021'};
sex=[2,2,2,2,2,2,2,2,2,2,2,2]; %2 male, 1 female
genotype=[1,1,1,1,3,3,3,3,2,2,2,2]; %1=85% WR, 2=CA, 3=95% WR
rawdata_folder='\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\Behavior\CitricAcid\bpodraw\';
GNG_protocol='LickingGNG_Normal_cam'; %'LickingGNG_Reverse_cam';
save_name='data_GNG_Normal'; %'data_GNG_Reverse';
save_folder=rawdata_folder;

Animals=[CA2]; %[APP10to12]; [CA]; [KF, AM, FZ ];

%______DO NOT MODIFY FROM HERE ON________________
clear data SessionData
fl=0;
for fi=1:length(Animals)
    fl=fl+1;
    folder=rawdata_folder; %'\\pbs-kk-labserv.win.ad.jhu.edu\Data18\Sharlen\Raw_Bpod_Data\';%'\\pbs-kk-labserv.win.ad.jhu.edu\Data5\LabData5\Sharlen\Bpod_raw\'; %'/Users/shar/Documents/APP_project/';
    foldertouse=[folder cell2mat(Animals(fl))  filesep GNG_protocol filesep]; %[folder cell2mat(Animals(fl))  filesep 'LickingGNG_Quarter' filesep]; % 'LickingGNG_Normal [folder cell2mat(Animals(fl)) '/LickingGNG_Quarter/'];
    if ~exist([foldertouse 'SessionsOut'])
        mkdir(foldertouse, 'SessionsOut') %makes directory to add here files that have less than 10 trials to be excluded from analysis
    end 
    cd([foldertouse 'Session Data' filesep]);
    matfiles2 = dir('*.mat'); %finds all GNG sessions for the selected animal
    
    %to make sure files with less than 10 trials are not taken into
    %consideration
    for maf=1:length(matfiles2)
        clear namemat
        namemat=matfiles2(maf).name;
        load(namemat);
        if SessionData.nTrials<10
            movefolder= ([foldertouse 'SessionsOut']);
            movefile(namemat, movefolder) %moves files to sessionOut folder per animal
        end
        cd([foldertouse 'Session Data' filesep]);
    end
    
    clear dA dA2 matfiles2
    matfiles2 = dir('*.mat'); %finds all GNG sessions for the selected animal
    [dA]=bpod_filesbinder(matfiles2); %finds index of files that were recorded on the same day to bind them
    
    %For GNG
    for fii= 1:length(matfiles2) %for all files from selected animal
        namemat=matfiles2(fii).name;
        load(namemat); %loads file by file to perform analysis below
        expday=str2num(namemat(end-18:end-11));
        Trials=SessionData.nTrials;
        Trials2use=1:SessionData.nTrials;
        
        for li=Trials2use 
            xvect=SessionData.RawEvents.Trial{1,Trials2use(li)};
            TrialStimOnset(li,:)=xvect.States.Stimulus(1);
            TrialStimDuration(li,:)=xvect.States.Stimulus(2)-TrialStimOnset(li,:);
            TrialToneGO(li,:)=[SessionData.TrialSettings(li).GUI.SinWaveFreqGo];
            TrialToneNOGO(li,:)=[SessionData.TrialSettings(li).GUI.SinWaveFreqNoGo];
            StimOnset(li,:)=SessionData.RawEvents.Trial{1,li}.States.Stimulus(1);
            RewardTime(li,:)=SessionData.RawEvents.Trial{1,li}.States.OpenValve(1);
            if isfield(xvect.Events, 'Port1In')==1 %if there was a lick
                Licks{li}=xvect.Events.Port1In; %absolute timing of all licks
                LickNum(li,:)=length(xvect.Events.Port1In);
                LickLickOnset{li}=xvect.Events.Port1In-xvect.Events.Port1In(1); %get lick timing aligned to the first lick
                TrialLickOnset{li}=xvect.Events.Port1In-TrialStimOnset(li,:); %get lick onset time (Drinking(1)-)
            else
                Licks{li}=NaN;
                LickNum(li,:)=NaN;
                LickLickOnset{li}=NaN;
                TrialLickOnset{li}=NaN; %if there was no lick its a NaN
            end
            if ~isnan(xvect.States.Punish(1)) %1=Hit, 2=Miss, 3=CR, 4=FA
                TrialResp(li,:)=4; %FA
            elseif ~isnan(xvect.States.CorrectReject(1))
                TrialResp(li,:)=3; %CR
            elseif ~isnan(xvect.States.Miss(1))
                TrialResp(li,:)=2; %Misses
            elseif ~isnan(xvect.States.Drinking(1))
                TrialResp(li,:)=1; %Hits
            end
            %to get time to first lick for the trials to use
            if isempty(find(TrialLickOnset{li}>0,1))
                FirstLickOnset(li,:)=NaN;
            else
                FirstLickOnset(li,:)=TrialLickOnset{li}(find(TrialLickOnset{li}>0,1));
            end
        end
        
        %For context
        ReinfTarget=find(SessionData.TrialTypes(Trials2use)==1); %finds Reinforcement target trials
        ReinfFoil=find(SessionData.TrialTypes(Trials2use)==2); % finds Reinforcement foil trials
        ProbeTarget=find(SessionData.TrialTypes(Trials2use)==3); % finds Probe target trials
        ProbeFoil=find(SessionData.TrialTypes(Trials2use)==4); %finds Probe foil trials
        ReinfAll=[ReinfTarget ReinfFoil];
        ProbeAll=[ProbeTarget ProbeFoil];
        
        %for Lick number
        TLicksRT=nansum(LickNum(ReinfTarget)); %total amount of licks during reinforcement target tone
        TLicksRF=nansum(LickNum(ReinfFoil)); %total amount of licks during reinforcement foil tone
        TLicksPT=nansum(LickNum(ProbeTarget)); %total amount of licks during probe target tone
        TLicksPF=nansum(LickNum(ProbeFoil)); %total amount of licks during probe foil tone
        
        %for Lick Onset
        LickOnset.RT=FirstLickOnset(ReinfTarget);
        LickOnset.RF=FirstLickOnset(ReinfFoil);
        LickOnset.PT=FirstLickOnset(ProbeTarget);
        LickOnset.PF=FirstLickOnset(ProbeFoil);
        
        %for Reinforcement context
        ReinfHit=length(intersect(ReinfAll, find(TrialResp==1)));
        ReinfMiss=length(intersect(ReinfAll, find(TrialResp==2)));
        ReinfCR=length(intersect(ReinfAll, find(TrialResp==3)));
        ReinfFA=length(intersect(ReinfAll, find(TrialResp==4)));
        
        %for Probe context
        ProbeHit=length(intersect(ProbeAll, find(TrialResp==1)));
        ProbeMiss=length(intersect(ProbeAll, find(TrialResp==2)));
        ProbeCR=length(intersect(ProbeAll, find(TrialResp==3)));
        ProbeFA=length(intersect(ProbeAll, find(TrialResp==4)));
        if (ProbeHit==0 & ProbeMiss==0 & ProbeCR==0 & ProbeFA==0)==1 %make NaNs in probe context if it was not present
            ProbeHit=NaN; ProbeMiss=NaN; ProbeCR=NaN; ProbeFA=NaN;
        end
        
        [dp_r, c_r]=dprime_criterion_calc(ReinfHit, ReinfMiss, ReinfCR, ReinfFA);
        [dp_p, c_p]=dprime_criterion_calc(ProbeHit, ProbeMiss, ProbeCR, ProbeFA);
        
        %to save data
        data(fl).mouse=cell2mat(Animals(fl));
        data(fl).expday{fii}=expday;
        data(fl).ntrials{fii}=length(Trials2use);
        data(fl).trialresp{fii}=TrialResp';
        data(fl).reinforced{fii}=[ReinfHit, ReinfMiss, ReinfCR, ReinfFA];
        data(fl).probe{fii}=[ProbeHit, ProbeMiss, ProbeCR, ProbeFA];
        data(fl).trialhist{fii}=SessionData.TrialTypes(Trials2use);
        data(fl).firstlickonset{fii}=LickOnset;
        data(fl).licksnum{fii}=[TLicksRT, TLicksRF, TLicksPT, TLicksPF]; 
        data(fl).licksall{fii}=Licks; %absolute timing of all licks
        data(fl).lickstime{fii}=TrialLickOnset; % all licks per trial aligned to the sound onset
        data(fl).lickstimefirst{fii}=LickLickOnset; %all licks per trial aligned to the first lick
        data(fl).rewardtime{fii}=RewardTime; %reward onset absolute timing
        data(fl).stimonset{fii}=StimOnset; %stimulus onset absolute timing
        data(fl).tonego{fii}=unique(TrialToneGO);
        data(fl).tonenogo{fii}=unique(TrialToneNOGO);
        data(fl).stimduration{fii}=TrialStimDuration;
        data(fl).dprime{fii}=[dp_r;dp_p];
        data(fl).criterion{fii}=[c_r;c_p];
        data(fl).sex=sex(fl);
        data(fl).genotype=genotype(fl);
        data(fl).dob=char(dob(fl));        
        
        clearvars -except data Animals fl folder dA dA2 matfiles2 dob sex genotype rawdata_folder GNG_protocol save_name save_folder
    end

    %Check more than 1 sesion on the same day (when Bpod crashed) and
    %change the data structure
    if iscell(dA) %~isnan(dA)
        for di=1:length(dA)
            predatt=[]; predatt2=[]; predatt3=[]; predatt4=[]; predatre=[]; predatpr=[]; predatln=[]; predatlonset=[]; predattrh=[]; predattrr=[]; predatlo=[]; predatlo2=[]; predatlo3=[]; predatsti=[]; predattonego=[]; predattonenogo=[];
            d=cell2mat(dA(di));
            for dii=1:length(d)
                predatt=[predatt cell2mat(data(fl).ntrials(d(dii)))]; %used trial number
                predatt3=[predatt3; cell2mat(data(fl).stimonset(d(dii)))]; %stim onset
                predatt4=[predatt4; cell2mat(data(fl).rewardtime(d(dii)))]; %reward onset
                predatre(dii,:)=cell2mat(data(fl).reinforced(d(dii))); %reinforced
                predatpr(dii,:)=cell2mat(data(fl).probe(d(dii))); %probe
                predatln=[predatln; cell2mat(data(fl).licksnum(d(dii)))]; %number of licks per condition 
                predatlonset=[predatlonset cell2mat(data(fl).firstlickonset(d(dii)))]; %lick onset time per condition
                predattrh=[predattrh cell2mat(data(fl).trialhist(d(dii)))]; %trial history
                predattrr=[predattrr cell2mat(data(fl).trialresp(d(dii)))]; %trial response history
                predattonego=[predattonego cell2mat(data(fl).tonego(d(dii)))];
                predattonenogo=[predattonenogo cell2mat(data(fl).tonenogo(d(dii)))];
                LO=data(fl).lickstime(d(dii));
                LO2=data(fl).lickstimefirst(d(dii));
                LO3=data(fl).licksall(d(dii));
                predatlo=horzcat(predatlo, LO{1,1}); %licks aligned to the sound
                predatlo2=horzcat(predatlo2, LO2{1,1}); %licks aligned to first lick
                predatlo3=horzcat(predatlo3, LO3{1,1});% all licks
                predatsti=[predatsti; cell2mat(data(fl).stimduration(d(dii)))]; %stimulus duration
            end
            
            if length(unique(predattonego))==1
                data(fl).tonego(d(1))={unique(predattonego)};
                data(fl).tonenogo(d(1))={unique(predattonenogo)};
            else
                disp(['Error in tone identity file ' matfiles2(dii).name])
                return
            end 
            
            data(fl).ntrials(d(1))={sum(predatt)}; %rewrite trial number in the min day
            data(fl).reinforced(d(1))={nansum(predatre)}; %rewrite reinforced number in the min day
            if sum(sum(isnan(predatpr)))==size(predatpr,1)*size(predatpr,2) %if all numbers in matrix are nans, we need to perform a sum
                data(fl).probe(d(1))={sum(predatpr)}; %rewrite reinforced number in the min day
            else %if there are numbers mixed with nans, we need to perform a nansum
                data(fl).probe(d(1))={nansum(predatpr)}; %rewrite probe number in the min day
            end
            
%             fnam=fieldnames(data(fl).firstlickonset{1,1}(1));
%             for fio=1:length(fnam)
%                 data(fl).firstlickonset{1,1}.(cell2mat(fnam(fio)))
%             end
            data(fl).trialhist(d(1))={predattrh}; %rewrite trial history in the min day
            data(fl).trialresp(d(1))={predattrr}; %rewrite trial response history in the min day
            data(fl).licksnum(d(1))={sum(predatln)};
            data(fl).lickstime(d(1))={predatlo}; %rewrite lick timing aligned to sound onset in the min day
            data(fl).lickstimefirst(d(1))={predatlo2}; %rewrite lick timing aligned to first lick in the min day
            data(fl).licksall(d(1))={predatlo3}; %rewrite lick timing of all licks unaligned
            data(fl).stimduration(d(1))={predatsti}; %rewrite stimulus duration in the min day
            data(fl).stimonset(d(1))={predatt3}; %rewrite stimulus onset in the min day
            data(fl).rewardtime(d(1))={predatt4}; %rewrite reward onset in the min day
            
            %redo dprime and criterion analysis
            [dp_r, c_r]=dprime_criterion_calc(cell2mat(data(fl).reinforced(d(1))));
            [dp_p, c_p]=dprime_criterion_calc(cell2mat(data(fl).probe(d(1))));
            data(fl).dprime(d(1))={[dp_r;dp_p]};
            data(fl).criterion(d(1))={[c_r;c_p]};
            clear predatt predatre predatpr predattrh predatlo predatlo2 predatlo3 predatsti dp_r dp_p c_r c_p LO LO2 LO3 predatlonset predatln predatt2 predatt3 predatt4 predattonego predattonenogo
        end
        
        %eliminate and shift the experiment days
        f=fields(data(fl));
        if length(dA)>=1 %in case there is only one experiment repetition for all training days
            for k=2:numel(f)-3 %for all variables except mouse number (1) sex, genotype and dob (20-22)
                for m=fliplr(1:length(dA)) %to start eliminating the not necessary experiments that have been incorporated, from the back to the front to avoid shifting the matrix
                    dA2=dA{m};
                    for mi=fliplr(2:size(dA2,1)) %to start eliminating the not necessary experiments that have been incorporated, from the back to the front to avoid shifting the matrix
                        data(fl).(f{k})(dA2(mi))=[]; %eliminates
                    end
                end
            end
        end
        clear d m LO k f
    end
end

% cd(folder)
% if ~exist([folder 'Analysis'], 'dir')
%    mkdir(folder, 'Analysis');
% end
% save([folder 'Analysis/data_GNGQuarter_' char(dayt) '.mat'], 'data')


%to save 
foldersave=save_folder;%'\\pbs-srv2.win.ad.jhu.edu\kishorelab\Sharlen\APP_Project\Behavior\'; %\KishoreLab\Sharlen\Behavior\CitricAcid\';
cd(foldersave)

if ~exist([foldersave 'Analysis'], 'dir')
   mkdir(foldersave, 'Analysis');
end

dayt=datetime('now', 'Format', 'yyyyMMdd''_''HHmmss');
save([foldersave 'Analysis' filesep save_name '_' char(dayt) '.mat'], 'data') %'data_GNGNormal_'

