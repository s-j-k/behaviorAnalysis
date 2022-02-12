%%%% Licking GNG Analysis %%%%%
%%%% Rashi Monga: 1 June 2021%%
%%%% Date of edits: 12 September 2021

%%% The following code has been adapted from Sarah's phasicAnalysis.m
%%% and has been modified to account for the Normal Licking GNG analysis

%%% Outcomes: a)Action Rate analysis of each individual animal i.e, Hit
%%% Rate and False Alarm Rate
%%% B) D' analysis for each individual animal
%%% C) Lick rate analysis
%%% D) Lick PSTH
%%% E) Lick latency
%%% F) Lick latency over the trials
%%% G) All the above analysis across animals as well
%%% If you are using BPOD1, 2, 3 then Rig Number is 5
%%% but if you are using BPOD4, then Rig Number is 6- this is all based on
%%% how ports are assigned for a given function in the BPOD

%%% Please note that this code picks up data from the excel file exactly
%%% from some columns so don't change the columns; add extra info to later
%%% columns

clear all
close all
% dataPath='T:\su\DATA\behaviorData\';
dataPath='S:\su\';
cd(dataPath)
% prompt1 = 'Enter the cohort';
cohortName = 'BPOD4';
cd (cohortName)
% prompt2 = 'Enter the file';
filename = 'dataframe';

cohortIn=cohortName; % find the animals (rows) that are part of the cohort you want
[num,txt,table] =xlsread((strcat(filename, '.xlsx')));
pp = strcmpi(cohortIn,table(:,13));
rowNum = find(pp==1);%Get Row number

animals = table(rowNum,1);
rig = table(2:end,7);
Code = table(2:end, 11);
ReinfTrials = table(2:end,12);
%virus=table(2:end,3);


for i= 1:length(animals)
    clear session context trial tone outcome licking date time sessionDay trialDay contextDay outcomeDay
    animalPath = [dataPath cohortName '\' cell2mat(animals(i)) '\' Code{i} '\Session Data\'];
    cd (animalPath);
    subj=cell2mat(animals(i,:));
    %     subjectPath = strcat('A:\sarah\Animals\',subj, '\*.mat'); % may have to alter file path here and on line 9
    file = dir(animalPath); %  specify files
    numFiles = length(file); % number of files
    %     filenames = strin; % create empty string array for file names
    Day = 1;
    for sess = 1:numFiles
        %         filenames(sess) = strcat('A:\sarah\Animals\',subj,'\', file(sess).name); % increment through each file
        filenamelength = length(file(end).name)-11;
        if ~isnan(strfind(file(sess).name, animals{i}))
            load(file(sess).name); % load each file
            day{sess}=datetime(SessionData.Info.SessionDate,'Format','yyyyMMdd');
            time{sess}=datetime(SessionData.Info.SessionStartTime_UTC, 'Format', 'HHmm');
            trial{sess} = num2cell(1:SessionData.nTrials);
            for tri = 1:SessionData.nTrials  % increment through all the trials in one given session
                session{sess}{tri} = sess; % SESSION NUMBER
                context{sess}{tri} = SessionData.TrialSettings(1).context(tri);
                if rig{i} ==5 && isfield(SessionData.RawEvents.Trial{1,tri}.Events,'Port1In')
                    licking{sess}{tri} = SessionData.RawEvents.Trial{1,tri}.Events.Port1In'-SessionData.RawEvents.Trial{1,tri}.States.Stimulus(1); % licks adjusted for stim
                elseif rig{i} ==5 && ~isfield(SessionData.RawEvents.Trial{1,tri}.Events,'Port1In')
                    licking{sess}{tri} = NaN; end
                if rig{i} ==6 && isfield(SessionData.RawEvents.Trial{1,tri}.Events,'Port2In')
                    licking{sess}{tri} = SessionData.RawEvents.Trial{1,tri}.Events.Port2In'-SessionData.RawEvents.Trial{1,tri}.States.Stimulus(1); % licks adjusted for stim
                elseif rig{i} ==6 && ~isfield(SessionData.RawEvents.Trial{1,tri}.Events,'Port2In')
                    licking{sess}{tri} = NaN; end
                % CONTEXT, 0=probe 1=reinf,2=catch targ,3=catch foil
                % celine's code: tril type: 1= targ, reinf; 2=foil,
                % reinf;3=targ, probe; 4=foil, probe; 5=targ,opto; 6=foil,opto
                % 7=targ, catch light off; 8=foil, catch light off, 9=targ,
                % catch light on; 10=foil, catchlight on.
                % OUTCOME (1=H,2=M,3=F,4=C) & LICK LATENCY & TONE (target=1,foil=2)
                if ~isnan(SessionData.RawEvents.Trial{1,tri}.States.OpenValve) % hit state
                    outcome{sess}{tri} = 1; tone{sess}{tri} = 4756;
                elseif ~isnan(SessionData.RawEvents.Trial{1,tri}.States.CorrectReject) %correct reject state
                    outcome{sess}{tri} = 4; tone{sess}{tri} = 8000;
                end
                if ~ismember(SessionData.TrialTypes(tri),[9 10 7 8])
                    if ~isnan(SessionData.RawEvents.Trial{1,tri}.States.Miss) %miss state, no lick
                        outcome{sess}{tri} = 2; tone{sess}{tri} = 4756;
                    end
                end
                if ismember(SessionData.TrialTypes(tri),[1 2 3 4 7 8])% || contains(subj,otherrig) % all but light on
                    if ~isnan(SessionData.RawEvents.Trial{1,tri}.States.Punish) %false alarm state
                        outcome{sess}{tri} = 3; tone{sess}{tri} = 8000;
                    end
                else % for light on, foil trials
                    if ~isnan(SessionData.RawEvents.Trial{1,tri}.States.PunishLightON) %false alarm state
                        outcome{sess}{tri} = 3; tone{sess}{tri} = 8000;
                    end
                end
            end %through all the trials in one session
            
            
            %                 sessionDay{1} = session{1,1};
            %     trialDay {1} = trial{1,1};
            %     contextDay {1} = context{1,1};
            %     outcomeDay {1} = outcome{1,1};
            %     lickDay {1} = licking{1,1};
            
            if sess==1 && ~isnan(strfind(file(sess).name, animals{i}))
                sessionDay{1} = cell2mat(session{1,1});
                trialDay {1} = cell2mat(trial{1,1});
                contextDay {1} = cell2mat(context{1,1});
                outcomeDay {1} = cell2mat(outcome{1,1});
                lickDay {1} = (licking{1,1});
            else
                
                if strncmpi(file(sess).name, file(sess-1).name, filenamelength)==1
                    sessionDay{Day} = cell2mat(horzcat(sessionDay{1, Day}, session{1, sess}));
                    trialDay {Day} = cell2mat(horzcat(trialDay{1,Day}, trial{1, sess}));
                    contextDay {Day} = cell2mat(horzcat(contextDay{1, Day}, context{1, sess}));
                    outcomeDay {Day} = cell2mat(horzcat(outcomeDay{1,Day}, outcome{1, sess}));
                    lickDay {Day} =(horzcat(lickDay{1,Day}, licking{1, sess}));
                else
                    Day = Day +1;
                    sessionDay{Day} = cell2mat(session{1, sess});
                    trialDay{Day} = cell2mat(trial{1, sess});
                    contextDay {Day}= cell2mat(context{1, sess});
                    outcomeDay {Day}=cell2mat(outcome{1, sess});
                    lickDay {Day}=(licking{1, sess});
                end
            end
        else
        end
    end
    datan(i).subj=animals{i};
    datan(i).gend=table{i+1,5};
    datan(i).strain=table{i+1,2};
    datan(i).virus=table{i+1,3};
    datan(i).exp=table{i+1,4};
    %     datan(i).session=horzcat(session{1,:});
    %     datan(i).trial=horzcat(trial{1,:});
    %      datan(i).context=horzcat(context{1,:});
    %     datan(i).tone=horzcat(tone{1,:});
    %     datan(i).outcome=horzcat(outcome{1,:});
    %     datan(i).lick=horzcat(licking{1,:});
    datan(i).sessDay = sessionDay;
    datan(i).trialDay = trialDay;
    datan(i).contextDay = contextDay;
    datan(i).outcomeDay = outcomeDay;
    datan(i).lickDay = lickDay;
    datan(i).date=horzcat(day{1,:});
    datan(i).time=horzcat(time{1,:});
    datan(i).trainer=table{i+1,6};
    datan(i).rig=table{i+1,7};
    %     datan(ani).age=table{ani+1,6};
    datan(i).use=table{i+1,10};
    %     datan(ani).note='NaN';
    disp(i);
    
    cd ../../..
end

for i = 1:length(datan)
    clear Trial Outcome Licks Context
    for m = 1:length(datan(i).sessDay)
        clear ExpContext ExpContextNaN ExpLick ExpOutcome ExpLick xlicks Exptrial
        if isempty(cell2mat((datan(i).sessDay(1,m))))==0
            if length(cell2mat((datan(i).sessDay(1,m))))~= (ReinfTrials{i}+20)
                ExpContextNaN = NaN.*(ones(1, ReinfTrials{i}+20));
                ExpContext = [ones(1, ReinfTrials{i}/2) zeros(1,20) ones(1, ReinfTrials{i}/2)];
                
                Exptrial =  NaN.*(ones(1, ReinfTrials{i}+20));
                Exptrial(1:length(cell2mat((datan(i).sessDay(1,m))))) = cell2mat(datan(i).trialDay(1,m));
                
                ExpOutcome = NaN.*(ones(1, ReinfTrials{i}+20));
                ExpOutcome(1:length(cell2mat((datan(i).sessDay(1,m))))) = cell2mat(datan(i).outcomeDay(1,m));
                
                ExpLick = cell(1, ReinfTrials{i}+20);
                xlicks = datan(i).lickDay(1,m);
                ExpLick (1:length(cell2mat((datan(i).sessDay(1,m)))))= xlicks{1, 1};
            else
                ExpContext = cell2mat(datan(i).contextDay(1,m));
                ExpOutcome = cell2mat(datan(i).outcomeDay(1,m));
                xlicks = datan(i).lickDay(1,m);
                ExpLick = xlicks{1,1};
                Exptrial = cell2mat(datan(i).trialDay(1,m));
            end
            Trial{1, m} = Exptrial;
            Outcome {1, m} = ExpOutcome;
            Licks{1, m} = ExpLick;
            Context{1, m} = ExpContext;           
        else
        end
        
    end
    datan(i).trials = num2cell(horzcat(Trial{1,:}));
    datan(i).outcome = num2cell(horzcat(Outcome{1,:}));
    datan(i).context = num2cell(horzcat(Context{1,:}));
    datan(i).lick =num2cell(horzcat(Licks{1,:}));
end

%%%% Day-wise analysis
for k = 1:length(datan)
    j = 1;
    for m = 1:length(datan(k).sessDay)
        clear Context Outcome Licks probeIndDay reinfIndDay reinfOutDay numHits numFA numHitsProbe numFAProbe probeIndDay probeOutDay probeOutTrial
        if isempty(cell2mat((datan(k).sessDay(1,m))))==0
            Context = cell2mat(datan(k).contextDay(1, m));
            Outcome = cell2mat(datan(k).outcomeDay(1,m));
            Licks = (datan(k).lickDay(1,m));
            %%% Day wise reinforced context
            reinfIndDay=find(Context==1);
            reinfOutDay = Outcome(reinfIndDay);
            numHits=length(find(reinfOutDay==1)); % find all hits
            numFA=length(find(reinfOutDay==3)); % find all fa
            DataDay(k).reinfhit(j)= numHits/(numHits+length(find(reinfOutDay==2)));
            DataDay(k).reinffa(j)= numFA/(numFA+length(find(reinfOutDay==4)));
            
            %%% Day wise probe context
            probeIndDay=find(Context==0); % find ind for probe
            probeOutDay=Outcome(probeIndDay);
            numHitsProbe = length(find(probeOutDay ==1));
            numFAProbe = length(find(probeOutDay==3));
            DataDay(k).probehit(j)= numHitsProbe/(numHitsProbe+length(find(probeOutDay==2)));
            DataDay(k).probefa(j)= numFAProbe/(numFAProbe+length(find(probeOutDay==4)));
            tstep2 = 20;
            dprimeRate_p(k).probehit=DataDay(k).probehit;
            dprimeRate_p(k).probefa=DataDay(k).probefa;
            dprimeRate_p(k).probehit(dprimeRate_p(k).probehit==1)=1-1/(2*tstep2); % normalize for dprime
            dprimeRate_p(k).probehit(dprimeRate_p(k).probehit==0)=1/(2*tstep2); % normalize for dprime
            dprimeRate_p(k).probefa(dprimeRate_p(k).probefa==1)=1-1/(2*tstep2); % normalize for dprime
            dprimeRate_p(k).probefa(dprimeRate_p(k).probefa==0)=1/(2*tstep2); % normalize for dprime
            DataDay(k).probemiss(j)=1-DataDay(k).probehit(j);
            DataDay(k).probecr(j)=1-DataDay(k).probefa(j);
            DataDay(k).dprimeprobe(j) = (norminv(full(dprimeRate_p(k).probehit(j)))) - (norminv(full(dprimeRate_p(k).probefa(j))));
            
            j = j+1;
        else
        end
    end
end


for subj = 1:length(datan) % go through each subj
    clear block bins_Reinf bins_Pr bins_licks bins_Pr reinfOut probeOut probeOutTrial
    probeInd=find(cellfun(@(x)isequal(x,0),datan(subj).context)==1); % find ind for probe
    reinfInd=find(cellfun(@(x)isequal(x,0),datan(subj).context)==0); % find ind for reinf
    %     reinfOut=datan(subj).outcome; reinfOut(probeInd)={NaN}; % change probe to NaN
    %     probeOut=datan(subj).outcome(probeInd);
    reinfOut=datan(subj).outcome(reinfInd);
    probeOut=datan(subj).outcome(probeInd);
    probeOutTrial = datan(subj).outcome;
    probeOutTrial(reinfInd)= {NaN};
    
    tstep=70; % blocks of 100 for the block design analysis
    s1=size(reinfOut,2);
    block=s1-mod(s1,tstep);
    bins_Reinf=reshape(reinfOut(1:block),tstep,[]); %bins of 100 for reinf
    ReinfLicks = datan(subj).lick(reinfInd);
    bins_licks = reshape(ReinfLicks(1:block), tstep, []);
    dataNew(subj).bins_Reinf = bins_Reinf;
    dataNew(subj).bins_licks = bins_licks;
    dataNew(subj).probeTrials = probeOutTrial;
    for j=1:size(bins_Reinf,2)
        numHits=sum(cellfun(@(x)isequal(x,1),bins_Reinf(:,j))); % find all hits
        numFA=sum(cellfun(@(x)isequal(x,3),bins_Reinf(:,j))); % find all fa
        dataNew(subj).reinfhit(j)= numHits/(numHits+sum(cellfun(@(x)isequal(x,2),bins_Reinf(:,j))));
        dataNew(subj).reinffa(j)= numFA/(numFA+sum(cellfun(@(x)isequal(x,4),bins_Reinf(:,j))));
        dprimeRate_r(subj).reinfhit=dataNew(subj).reinfhit;
        dprimeRate_r(subj).reinffa=dataNew(subj).reinffa;
        dprimeRate_r(subj).reinfhit(dprimeRate_r(subj).reinfhit==1)=1-1/(2*tstep); % normalize for dprime
        dprimeRate_r(subj).reinfhit(dprimeRate_r(subj).reinfhit==0)=1/(2*tstep); % normalize for dprime
        dprimeRate_r(subj).reinffa(dprimeRate_r(subj).reinffa==1)=1-1/(2*tstep); % normalize for dprime
        dprimeRate_r(subj).reinffa(dprimeRate_r(subj).reinffa==0)=1/(2*tstep); % normalize for dprime
        dataNew(subj).reinfmiss(j)=1-dataNew(subj).reinfhit(j);
        dataNew(subj).reinfcr(j)=1-dataNew(subj).reinffa(j);
        dataNew(subj).dprimereinf(j) = (norminv(full(dprimeRate_r(subj).reinfhit(j)))) - (norminv(full(dprimeRate_r(subj).reinffa(j))));
        %         dataNew(subj).HitLicks(each)= datan(subj).lick(cellfun(@(x)isequal(x,1),bins(:,each)));
        idxHit = find((cellfun(@(x)isequal(x,1),bins_Reinf(:,j)))==0);
        dataNew(subj).HitLicks{j}= bins_licks(:, j);
        dataNew(subj).HitLicks{j}(idxHit)= {NaN};
        
        idxFA = find((cellfun(@(x)isequal(x,3),bins_Reinf(:,j)))==0);
        dataNew(subj).FALicks{j}= bins_licks(:, j);
        dataNew(subj).FALicks{j}(idxFA)= {NaN};
        
        idxCR = find((cellfun(@(x)isequal(x,4),bins_Reinf(:,j)))==0);
        dataNew(subj).CRLicks{j}= bins_licks(:,j);
        dataNew(subj).CRLicks{j}(idxCR)= {NaN};
    end
    %     tstep2 = 20;
    %     s2 = size(probeOutTrial, 2);
    %     block2=s2-mod(s2,tstep2);
    %     bins_Pr=reshape(probeOutTrial(1:block2),tstep2,[]); % bins of 100 for probe
    %     dataNew(subj).bins_Pr = bins_Pr;
    %     for j=1:size(bins_Pr,2)
    %         numHits=sum(cellfun(@(x)isequal(x,1),bins_Pr(:,j))); % find all hits
    %         numFA=sum(cellfun(@(x)isequal(x,3),bins_Pr(:,j))); % find all fa
    %         dataNew(subj).probehit(j)= numHits/(numHits+sum(cellfun(@(x)isequal(x,2),bins_Pr(:,j))));
    %         dataNew(subj).probefa(j)= numFA/(numFA+sum(cellfun(@(x)isequal(x,4),bins_Pr(:,j))));
    %         dprimeRate_p(subj).probehit=dataNew(subj).probehit;
    %         dprimeRate_p(subj).probefa=dataNew(subj).probefa;
    %         dprimeRate_p(subj).probehit(dprimeRate_p(subj).probehit==1)=1-1/(2*tstep2); % normalize for dprime
    %         dprimeRate_p(subj).probehit(dprimeRate_p(subj).probehit==0)=1/(2*tstep2); % normalize for dprime
    %         dprimeRate_p(subj).probefa(dprimeRate_p(subj).probefa==1)=1-1/(2*tstep2); % normalize for dprime
    %         dprimeRate_p(subj).probefa(dprimeRate_p(subj).probefa==0)=1/(2*tstep2); % normalize for dprime
    %         dataNew(subj).probemiss(j)=1-dataNew(subj).probehit(j);
    %         dataNew(subj).probecr(j)=1-dataNew(subj).probefa(j);
    %         dataNew(subj).dprimeprobe(j) = (norminv(full(dprimeRate_p(subj).probehit(j)))) - (norminv(full(dprimeRate_p(subj).probefa(j))));
    %     end
    % FIX MAX PROBE
    %     dataNew(subj).dprimeprobe(find(dataNew(subj).dprimeprobe==max(dataNew(subj).dprimeprobe)):length(dataNew(subj).dprimeprobe)) = max(dataNew(subj).dprimeprobe);
end



SingleAnimalPlots3(datan, dataNew, animals, DataDay);
% AcrossAnimalPlots(datan, dataNew, animals);

