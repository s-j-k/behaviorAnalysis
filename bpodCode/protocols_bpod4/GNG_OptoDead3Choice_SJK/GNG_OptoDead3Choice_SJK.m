function GNG_OptoDead3Choice_SJK

% Parameters that can be suitable to change are lines 27 and 256.
% currently coded for opto on 10% of trials
global BpodSystem
global Amount

%% ENTER HERE THE SESSION NUMBER OF THE MOUSE %%
trainingday = 10;

%% Create trial manager object 
TrialManager = TrialManagerObject;

%% Define parameters
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    S.GUI.RewardAmount = 3.0; % ul % due to calibration for output
    S.GUI.SoundDuration = 0.1; % duration of sound
%     S.GUI.SinWaveFreqGo = 8000; % Frequency of go cue
%     S.GUI.SinWaveFreqNoGo = 4756; % Frequency of no-go %ue
    S.GUI.SinWaveFreqGo =  4756;% Frequency of go cue
    S.GUI.SinWaveFreqNoGo = 8000; % Frequency of no-go cue    
%     
    S.GUIPanels.Sound = {'SinWaveFreqGo', 'SinWaveFreqNoGo', 'SoundDuration'}; % Labels for sound panel
end
Amount = S.GUI.RewardAmount;%actual calibration amount (ul from g)

%% Go and NoGo factors for sound amplitudes

% Su calibration, target 70dB, 2/15/2024
% go_factor = 0.6; % if go = 8000, dB = 66-67dB %0.6  
% nogo_factor = 1; % if nogo = 4756, dB = 65-66dB %1
go_factor = 1; % if go = 4756, dB = 66-67dB %1
nogo_factor = 0.5; % if nogo = 8000, dB = 65-66dB %0.6

%% Define trials

% Define default values here
MaxTrials = 400; % max trials (for display purposes)
ntrials = 300;
reinf = [1 MaxTrials]; 
% probe1 = [301 302]; % Input trials for first probe block (20 trials min)
probe1 = [141 160];
opto = 0.5; % proportion of reinforced trials
% opto = 0.4;
% catchh = 0.3; % proportion of reinforced trials. If opto>0, opto*light-on, (1-opto)*light-off
% opto= 0;
catchh=0;
onlyoptoincatch = true;
fiftyfiftyincatch = false;

onlyoptoinseconddead = true; % this will replace 1/2 of the full trial inactivation 
onlyoptoinfirstdead = false; % this will replace 1/2 of the full trial inactivation. 
% if one of these values is set to false, the optogenetic inactivation will
% occur in the middle of the dead period and last 100ms
% second dead happens 100ms into the dead period until the first lick
% first dead period occurs at the beginning of the dead period and lasts
% for 100ms

% Other default values
nToneCounterbalanced = 20; % 20 minimum
val = MaxTrials/nToneCounterbalanced;
nSwitchTargetTones = 2; % nb of target tones play after a switch (i.e. from Reinforced to Probe or vice versa)
nProbe = diff(probe1)+1;
if ~isempty(probe1)
    nReinf = diff(reinf)+1 - nProbe;
else
    nReinf = diff(reinf)+1;
end
target = 1;
foil = 2;  

% Despite having separate trial types still require different contexts on top of that because separate state machine needs to run with licktube out (multiple valves cannot be open at once)
% Generate context vector for the session (reinforced, probe, reinforced opto)
S.context = zeros(MaxTrials,1);
if ~isempty(reinf)
    S.context(reinf(1):reinf(2)) = 2;% 2 = reinf context, licktube in 
end
if ~isempty(probe1)
    S.context(probe1(1):probe1(2)) = 0;% 0 = probe context, licktube out
end

% Generate trial type vector for the session
nDaysFiveTarget = 5; % nb of trainingday with 5 target tones at the beginning of session
if trainingday <= nDaysFiveTarget
    nFirstTargetTones = 5;
    nToneCounterbalanced_firstbloc = 20; % nb of counterbalanced tones in the first bloc when 5 1st tones are TARGET 
else 
    nFirstTargetTones = 2;
    nToneCounterbalanced_firstbloc = 20; % nb of counterbalanced tones in the first bloc when 2 1st tones are TARGET 
end

% Generate initial vector:
initialvector = [foil*ones(nToneCounterbalanced/2,1);target*ones(nToneCounterbalanced/2,1)];
TrialTypes = zeros(nToneCounterbalanced,val);
for x = 1:val
    TrialTypes(:,x) = initialvector(randperm(nToneCounterbalanced));
end
TrialTypes = reshape(TrialTypes,MaxTrials,1);

% Now apply restrictions:
% 1) first 5 (or 2) tones are target (i.e counterbalanced the first 20-trial (or 10-trial) bloc accordingly)
tone_vector = [target*ones(nToneCounterbalanced_firstbloc/2-nFirstTargetTones,1);foil*ones(nToneCounterbalanced_firstbloc/2,1)];
first_bloc_tones = [target*ones(nFirstTargetTones,1);tone_vector(randperm(nToneCounterbalanced_firstbloc-nFirstTargetTones))];
TrialTypes(1:nToneCounterbalanced_firstbloc)= first_bloc_tones;

if ~isempty(probe1)
    % 2) Make sure that the first two tones after a switch are target tones
    tone_vector_swich = [target*ones(nToneCounterbalanced/2-nSwitchTargetTones,1);foil*ones(nToneCounterbalanced/2,1)];
    first_bloc_probe_tones = [target*ones(nSwitchTargetTones,1);tone_vector_swich(randperm(nToneCounterbalanced-nSwitchTargetTones))];
    TrialTypes(probe1(1):probe1(1)+nToneCounterbalanced-1) = first_bloc_probe_tones; 
    swich_active_tones = [target*ones(nSwitchTargetTones,1);tone_vector_swich(randperm(nToneCounterbalanced-nSwitchTargetTones))];
    TrialTypes(probe1(2)+1:probe1(2)+1+nToneCounterbalanced-1) = swich_active_tones; 
end

% If asked, select pseudo-random opto trials, every nToneCounterbalanced trials
% if opto no catch
if opto~=0 
    if ~isempty(probe1)
        blocReinf1 = [reinf(1) probe1(1)-1];
        nBlocReinf1 = diff(blocReinf1)+1;
        nOptoTrials1 = opto*nBlocReinf1;
        nOptoTrials1Conterbalenced = nOptoTrials1/(nBlocReinf1/nToneCounterbalanced);
        
        vec = TrialTypes(blocReinf1(1):blocReinf1(end));
        vec1 = find(vec==1);
        [~,w] = sort(randn(nToneCounterbalanced/2,(nBlocReinf1/nToneCounterbalanced)));
        idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
        idx = repelem(idx,nToneCounterbalanced/2,1);
        w = w+idx; 
        select = w(1:nOptoTrials1Conterbalenced/2,:);
        select = select(:);
        vec1octo = vec1(select);
        vec1octo = vec1octo(:);        
        vec2 = find(vec==2);
        [~,w] = sort(randn(nToneCounterbalanced/2,(nBlocReinf1/nToneCounterbalanced)));
        w = w+idx; 
        select = w(1:nOptoTrials1Conterbalenced/2,:);
        select = select(:);
        vec2octo = vec2(select);
        vec2octo = vec2octo(:);
        
        vecopto = sort([vec1octo;vec2octo]);
        S.context(vecopto) = 1;% 1 = opto context, licktube in 
                
        blocReinf2 = [probe1(2)+1 reinf(2)];
        nBlocReinf2 = diff(blocReinf2)+1;
        nOptoTrials2 = opto*nBlocReinf2;
        nOptoTrials1Conterbalenced = nOptoTrials2/(nBlocReinf2/nToneCounterbalanced);
                
        vec = TrialTypes(blocReinf2(1):blocReinf2(end));
        vec1 = find(vec==1);
        [~,w] = sort(randn(nToneCounterbalanced/2,(nBlocReinf2/nToneCounterbalanced)));
        idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
        idx = repelem(idx,nToneCounterbalanced/2,1);
        w = w+idx; 
        select = w(1:nOptoTrials1Conterbalenced/2,:);
        select = select(:);
        vec1octo = vec1(select);
        vec1octo = vec1octo(:);    
        vec1octo = vec1octo+probe1(2);
        vec2 = find(vec==2);
        [~,w] = sort(randn(nToneCounterbalanced/2,(nBlocReinf2/nToneCounterbalanced)));
        w = w+idx; 
        select = w(1:nOptoTrials1Conterbalenced/2,:);
        select = select(:);
        vec2octo = vec2(select);
        vec2octo = vec2octo(:);
        vec2octo = vec2octo+probe1(2);
        
        vecopto = sort([vec1octo;vec2octo]);
        S.context(vecopto) = 1;% 1 = opto context, licktube in
            
    else   % if there's no probe         
        nOptoTrials = opto*nReinf;
        nOptoTrials1Conterbalenced = nOptoTrials/(nReinf/nToneCounterbalanced);
        
        vec = TrialTypes;
        vec1 = find(vec==1);
        [~,w] = sort(randn(nToneCounterbalanced/2,(nReinf/nToneCounterbalanced)));
        idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
        idx = repelem(idx,nToneCounterbalanced/2,1);
        w = w+idx; 
        select = w(1:nOptoTrials1Conterbalenced/2,:);
        select = select(:);
        vec1octo = vec1(select);
        vec1octo = vec1octo(:);        
        vec2 = find(vec==2);
        [~,w] = sort(randn(nToneCounterbalanced/2,(nReinf/nToneCounterbalanced)));
        w = w+idx; 
        select = w(1:nOptoTrials1Conterbalenced/2,:);
        select = select(:);
        vec2octo = vec2(select);
        vec2octo = vec2octo(:);
        
        vecopto = sort([vec1octo;vec2octo]);
        S.context(vecopto) = 1;% 1 = opto context, licktube in 
    end
end

%-----CONVERT TrialTypes FROM 1/2 TO 1/2(REINF), 3/4(PROBE), 5/6(OPTO)-----
% CONTEXT (S.context) --> 0=PROBE 1=OPTO 2=NORMAL
% TRIALTYPES --> opto(G,NG) = 5,6; reinf(G,NG) = 1,2; probe(G,NG) = 3,4;
TrialTypes(S.context==0) = TrialTypes(S.context==0)+2; % probe context (=0)
TrialTypes(S.context==1) = TrialTypes(S.context==1)+4; % opto context (=1)

temp_TrialTypes = TrialTypes(1:ntrials);
%------Full Trial Opto Context to Tone-only Context or Choice context------
% sjk 2/8/2024
% of the opto trials, change some of them to a different trial type if
% you want tone-only or choice-only inactivation
% CONTEXT (S.context) --> 0=PROBE 1=OPTO 2=NORMAL
if onlyoptoinseconddead && onlyoptoinfirstdead %opto for Tone AND choice, no full trial
%     contextIdx=find(temp_S.context==1); % find context for opto  % you
%     don't need it (CD)
    goTypeIdx = find(temp_TrialTypes==5); % find opto Go tone
    nogoTypeIdx = find(temp_TrialTypes==6); % find opto No go tone
    % halve the counterbalanced trial vectors, first for go tones
    randomGoTypeIdx=randperm(length(goTypeIdx));
    firstHalfGoIdx = goTypeIdx(randomGoTypeIdx(1:end/2));
    [secondHalfGoIdx, secondHalfGoVal]= (setdiff(goTypeIdx, firstHalfGoIdx));
    TrialTypes(firstHalfGoIdx) = TrialTypes(firstHalfGoIdx)+6; % Trial type 11 Tone opto go
    TrialTypes(secondHalfGoIdx) = TrialTypes(secondHalfGoIdx)+8; % Trial type 13 Choice opto go
    % now for the no go tones
    randomNoGoTypeIdx=randperm(length(nogoTypeIdx));
    firstHalfNoGoIdx = sort(nogoTypeIdx(randomNoGoTypeIdx(1:end/2)));
    [secondHalfNoGoIdx, secondHalfNoGoVal]= (setdiff(nogoTypeIdx, firstHalfNoGoIdx));
    TrialTypes(firstHalfNoGoIdx) = TrialTypes(firstHalfNoGoIdx)+6; % Trial type 12 Tone opto go
    TrialTypes(secondHalfNoGoIdx) = TrialTypes(secondHalfNoGoIdx)+8; % Trial type 14 Choice opto go
    % update the context variables by getting the index of the trialtypes
    S.context(firstHalfGoIdx) = 5;% trialtype 11 and 12 = context 5 = tone opto
    S.context(secondHalfGoIdx) = 6;
    S.context(firstHalfNoGoIdx) = 5; % trialtype 13 and 14 = context 6 = choice opto
    S.context(secondHalfNoGoIdx) = 6;
    
elseif onlyoptoinseconddead && onlyoptoinfirstdead~=1 %opto for Tone AND Full Trial
    contextIdx=find(S.context==1); % find context for opto 
    % replace one counterbalanced half with onlyoptointone context 
    gotypeIdx = find(TrialTypes==5); % Go tone
    nogotypeIdx = find(TrialTypes==6); % no go
    goTypeIdx = find(temp_TrialTypes==5); % find opto Go tone
    nogoTypeIdx = find(temp_TrialTypes==6); % find opto No go tone
    % halve the counterbalanced trial vectors, first for go tones
    randomGoTypeIdx=randperm(length(goTypeIdx));
    firstHalfGoIdx = goTypeIdx(randomGoTypeIdx(1:end/2));
    TrialTypes(firstHalfGoIdx) = TrialTypes(firstHalfGoIdx)+6; % Trial type 11 Tone opto go
    randomNoGoTypeIdx=randperm(length(nogoTypeIdx));
    firstHalfNoGoIdx = sort(nogoTypeIdx(randomNoGoTypeIdx(1:end/2)));
    TrialTypes(firstHalfNoGoIdx) = TrialTypes(firstHalfNoGoIdx)+6; % Trial type 12 Tone opto go
    % update the context variables by getting the index of the trialtypes
    S.context(firstHalfGoIdx) = 5;% trialtype 11 and 12 = context 5 = tone opto
    S.context(firstHalfNoGoIdx) = 5;

    
elseif onlyoptoinfirstdead && onlyoptoinseconddead~=1   %opto for Choice AND Full Trial
    % replace the other counterbalanced half with onlyoptoinchoice context
        goTypeIdx = find(temp_TrialTypes==5); % find opto Go tone
    nogoTypeIdx = find(temp_TrialTypes==6); % find opto No go tone
    % halve the counterbalanced trial vectors, first for go tones
    randomGoTypeIdx=randperm(length(goTypeIdx));
    firstHalfGoIdx = goTypeIdx(randomGoTypeIdx(1:end/2));
    TrialTypes(firstHalfGoIdx) = TrialTypes(firstHalfGoIdx)+8; % Trial type 13 Choice opto go
    randomNoGoTypeIdx=randperm(length(nogoTypeIdx));
    firstHalfNoGoIdx = sort(nogoTypeIdx(randomNoGoTypeIdx(1:end/2)));
    TrialTypes(firstHalfNoGoIdx) = TrialTypes(firstHalfNoGoIdx)+8; % Trial type 14 Choice opto go
    % update the context variables by getting the index of the trialtypes
    S.context(firstHalfGoIdx) = 6;
    S.context(firstHalfNoGoIdx) = 6; % trialtype 13 and 14 = context 6 = choice opto
else
    % all full trial opto inactivation (: 
end


% Add 4th (and 5th, if opto>0) context: CATCH trials (context==3: reinf light-off; context==4: reinf light-on)
if opto~=0 && catchh~=0 && ~onlyoptoincatch && ~fiftyfiftyincatch
    % IF OPTO, CATCH, NO OPTO IN CATCH
    nReinfLightOff = sum(TrialTypes==1 | TrialTypes==2);
    ncatchLightOff = catchh*nReinfLightOff/nToneCounterbalanced;
    ncatchLightOff = floor(ncatchLightOff);
    
    vec1 = find(TrialTypes==1);
    if mod(ncatchLightOff,2)~=0, ncatchLightOff = ncatchLightOff-1; end
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    select = w(1:ncatchLightOff/2,:);
    select = select(:);    
    vec1catch = vec1(select);
    vec1catchLightOff = vec1catch(:);    
        
    vec2 = find(TrialTypes==2);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatchLightOff/2,:);
    select = select(:);    
    vec2catch = vec2(select);
    vec2catchLightOff = vec2catch(:);  
    
    TrialTypes(vec1catchLightOff) = 7;
    TrialTypes(vec2catchLightOff) = 8;
    
    veccatchLightOff = sort([vec1catchLightOff;vec2catchLightOff]);
    S.context(veccatchLightOff) = 3;
       
    % Light ON
    nReinfLightOn = sum(TrialTypes==5 | TrialTypes==6);
    ncatchLightOn = catchh*nReinfLightOn/nToneCounterbalanced;
    if ncatchLightOn<1, ncatchLightOn=1;
    else, ncatchLightOn = floor(ncatchLightOn); end
    
    vec1 = find(TrialTypes==5);    
    if mod(ncatchLightOn,2)~=0, ncatchLightOn = ncatchLightOn-1; end
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOn/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    select = w(1:ncatchLightOn/2,:);
    select = select(:);    
    vec1catch = vec1(select);
    vec1catchLightOn = vec1catch(:);    
        
    vec2 = find(TrialTypes==6);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOn/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatchLightOn/2,:);
    select = select(:);    
    vec2catch = vec2(select);
    vec2catchLightOn = vec2catch(:);  
    
    TrialTypes(vec1catchLightOn) = 9;
    TrialTypes(vec2catchLightOn) = 10;
    
    veccatchLightOn = sort([vec1catchLightOn;vec2catchLightOn]);
    S.context(veccatchLightOn) = 4;
elseif opto~=0 && catchh==0 && ~onlyoptoincatch && fiftyfiftyincatch
    % this one is going to be a problem... 2/8/24
    % IF OPTO, CATCH, 50/50 OPTO in CATCH
    nReinfLightOff = sum(TrialTypes==1 | TrialTypes==2);
    ncatchLightOff = catchh*nReinfLightOff/nToneCounterbalanced;
    ncatchLightOff = floor(ncatchLightOff);
    
    vec1 = find(TrialTypes==1);
    if mod(ncatchLightOff,2)~=0, ncatchLightOff = ncatchLightOff-1; end
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    select = w(1:ncatchLightOff/2,:);
    select = select(:);    
    vec1catch = vec1(select);
    vec1catchLightOff = vec1catch(:);    
        
    vec2 = find(TrialTypes==2);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatchLightOff/2,:);
    select = select(:);    
    vec2catch = vec2(select);
    vec2catchLightOff = vec2catch(:);  
    
    TrialTypes(vec1catchLightOff) = 7;
    TrialTypes(vec2catchLightOff) = 8;
    
    veccatchLightOff = sort([vec1catchLightOff;vec2catchLightOff]);
    S.context(veccatchLightOff) = 3;
       
    % Light ON (computed from light OFF)
    nReinfLightOn = nReinfLightOff;
    ncatchLightOn = catchh*nReinfLightOn/nToneCounterbalanced;
    if ncatchLightOn<1, ncatchLightOn=1;
    else, ncatchLightOn = floor(ncatchLightOn); end
    
    vec1 = find(TrialTypes==1);    
    if mod(ncatchLightOn,2)~=0, ncatchLightOn = ncatchLightOn-1; end
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(sum(TrialTypes==1 | TrialTypes==2)/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    select = w(1:ncatchLightOn/2,:);
    select = select(:);    
    vec1catch = vec1(select);
    vec1catchLightOn = vec1catch(:);    
        
    vec2 = find(TrialTypes==2);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(sum(TrialTypes==1 | TrialTypes==2)/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatchLightOn/2,:);
    select = select(:);    
    vec2catch = vec2(select);
    vec2catchLightOn = vec2catch(:);  
    
    TrialTypes(vec1catchLightOn) = 9;
    TrialTypes(vec2catchLightOn) = 10;
    
    veccatchLightOn = sort([vec1catchLightOn;vec2catchLightOn]);
    S.context(veccatchLightOn) = 4;
elseif opto~=0 && catchh~=0 && onlyoptoincatch && ~fiftyfiftyincatch
    
    % IF OPTO, CATCH, 100 OPTO IN CATCH 
    nReinfLightOff = sum(TrialTypes==1 | TrialTypes==2);
    ncatchLightOff = catchh*nReinfLightOff/nToneCounterbalanced;
    ncatchLightOff = floor(ncatchLightOff);
    
    vec1 = find(TrialTypes==1);
    if mod(ncatchLightOff,2)~=0, ncatchLightOff = ncatchLightOff-1; end
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    select = w(1:ncatchLightOff/2,:);
    select = select(:);    
    vec1catch = vec1(select);
    vec1catchLightOff = vec1catch(:);    
        
    vec2 = find(TrialTypes==2);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatchLightOff/2,:);
    select = select(:);    
    vec2catch = vec2(select);
    vec2catchLightOff = vec2catch(:);  
    
    TrialTypes(vec1catchLightOff) = 9;
    TrialTypes(vec2catchLightOff) = 10;
    
    veccatchLightOff = sort([vec1catchLightOff;vec2catchLightOff]);
    S.context(veccatchLightOff) = 4; % Ligh ON
elseif opto==0 && catchh~=0 && onlyoptoincatch &&  ~fiftyfiftyincatch
    
    % IF NO OPTO, CATCH, 100 OPTO CATCH
    ncatch = catchh*nReinf/nToneCounterbalanced;
    ncatch = floor(ncatch);
    if mod(ncatch,2)~=0, ncatch = ncatch-1; end
    
    vec1 = find(TrialTypes==1);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinf/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    select = w(1:ncatch/2,:);
    select = select(:);
    vec1catch = vec1(select);
    vec1catch = vec1catch(:);        
    vec2 = find(TrialTypes==2);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinf/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatch/2,:);
    select = select(:);
    vec2catch = vec2(select);
    vec2catch = vec2catch(:);

    TrialTypes(vec1catch) = 9;
    TrialTypes(vec2catch) = 10;
    
    veccatch = sort([vec1catch;vec2catch]);
    S.context(veccatch) = 4; % Ligh ON
elseif opto==0 && catchh~=0 && ~onlyoptoincatch && ~fiftyfiftyincatch
    
    % IF NO OPTO, CATCH, NO OPTO IN CATCH
    ncatch = catchh*nReinf/nToneCounterbalanced;
    ncatch = floor(ncatch);
    if mod(ncatch,2)~=0, ncatch = ncatch-1; end
    
    vec1 = find(TrialTypes==1);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinf/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    select = w(1:ncatch/2,:);
    select = select(:);
    vec1catch = vec1(select);
    vec1catch = vec1catch(:);        
    vec2 = find(TrialTypes==2);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinf/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatch/2,:);
    select = select(:);
    vec2catch = vec2(select);
    vec2catch = vec2catch(:);

    TrialTypes(vec1catch) = 7;
    TrialTypes(vec2catch) = 8;
    
    veccatch = sort([vec1catch;vec2catch]);
    S.context(veccatch) = 3; % Light OFF
    
    
    
%-----CHANGE BY ZZAW042912-----
% Designed for the first few days of 50/50 after acquisition; no opto trials paired with tone in behavior, only opto trials in catch trials

% CONTEXT (S.context) --> 0=PROBE 1=OPTO 2=NORMAL
% TRIALTYPES --> opto(G,NG) = 5,6; reinf(G,NG) = 1,2; probe(G,NG) = 3,4;
elseif opto==0 && catchh~=0 && ~onlyoptoincatch && fiftyfiftyincatch
    % IF NO OPTO, CATCH, 50/50 OPTO IN CATCH
    nReinfLightOff = sum(TrialTypes==1 | TrialTypes==2);
    ncatchLightOff = catchh/2*nToneCounterbalanced;
    ncatchLightOff = floor(ncatchLightOff);
    if mod(ncatchLightOff,2)~=0, ncatchLightOff = ncatchLightOff-1; end
    
    vec1 = find(TrialTypes==1);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    idx = 0:nToneCounterbalanced/2:(size(w,2)-1)*(nToneCounterbalanced/2);
    idx = repelem(idx,nToneCounterbalanced/2,1);
    w = w+idx; 
    
    % select light off catch target trial
    selectOff = w(1:ncatchLightOff/2,:);
    selectOff = selectOff(:);    
    vec1catch = vec1(selectOff);
    vec1catchLightOff = vec1catch(:);
    TrialTypes(vec1catchLightOff) = 7;
    S.context(vec1catchLightOff) = 3;
    
    % select light on catch target trial
    selectOn = w(ncatchLightOff/2+1:ncatchLightOff,:);
    selectOn = selectOn(:);    
    vec1catch = vec1(selectOn);
    vec1catchLightOn = vec1catch(:);
    TrialTypes(vec1catchLightOn) = 9;
    S.context(vec1catchLightOn) = 4;
    
    % SELECT FOIL TRIALS    
    vec2 = find(TrialTypes==2);
    [~,w] = sort(randn(nToneCounterbalanced/2,floor(nReinfLightOff/nToneCounterbalanced)));
    w = w+idx; 
    select = w(1:ncatchLightOff/2,:);
    select = select(:);    
    vec2catchOff = vec2(select);
    vec2catchLightOff = vec2catchOff(:); 
    TrialTypes(vec2catchLightOff) = 8;
    S.context(vec2catchLightOff) = 3;
    
    select = w(ncatchLightOff/2+1:ncatchLightOff,:);
    select = select(:);    
    vec2catchOn = vec2(select);
    vec2catchLightOn = vec2catchOn(:);  
    TrialTypes(vec2catchLightOn) = 10;
    S.context(vec2catchLightOn) = 4;
    
%------------------------------


end

% % TO CHECK the code
% Accumulate(TrialTypes)% 
% [sum(S.context==1) sum(S.context==2) sum(S.context==3) sum(S.context==4)]'

ctxx = S.context(1:ntrials);
trialT = TrialTypes(1:ntrials);

disp([num2str(sum(ctxx==2)) ' reinf, ' num2str(sum(trialT==1)) ' target, ' num2str(sum(trialT==2)) ' foil.']);
disp([num2str(sum(ctxx==0)) ' probe, ' num2str(sum(trialT==3)) ' target, ' num2str(sum(trialT==4)) ' foil.']);
disp([num2str(sum(ctxx==1)) ' Third Dead opto, ' num2str(sum(trialT==5)) ' target, ' num2str(sum(trialT==6)) ' foil.']);
disp([num2str(sum(ctxx==5)) ' Second Dead opto, ' num2str(sum(trialT==11)) ' target, ' num2str(sum(trialT==12)) ' foil.']);
disp([num2str(sum(ctxx==6)) ' First Dead Only opto, ' num2str(sum(trialT==13)) ' target, ' num2str(sum(trialT==14)) ' foil.']);
disp([num2str(sum(ctxx==3)) ' catch light-OFF, ' num2str(sum(trialT==7)) ' no target, ' num2str(sum(trialT==8)) ' no foil.']);
disp([num2str(sum(ctxx==4)) ' catch light-ON, ' num2str(sum(trialT==9)) ' no target, ' num2str(sum(trialT==10)) ' no foil.']);
% sum(TrialTypes==3 | TrialTypes==4)
% sum(TrialTypes==5 | TrialTypes==6)
% sum(TrialTypes==7 | TrialTypes==8)
% sum(TrialTypes==9 | TrialTypes==10)

% figure; plot(TrialTypes,'.');
% ylim([0 11]);
%%
BpodSystem.Data.TrialTypes = []; % The trial type of each trial completed will be added here.
tic; % starts timer
%% Initialize plots
BpodSystem.ProtocolFigures.OutcomePlotFig = figure('Position', [10 750 1900 400],'name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off'); % Initializes figure for Outcome plot
BpodSystem.GUIHandles.OutcomePlot = axes('Position', [.075 .3 .9 .6]); % Initializes axes for Outcome plot
TrialTypeOutcomePlot_Opto(BpodSystem.GUIHandles.OutcomePlot,'init',TrialTypes') % 'ntrials',MaxTrials); % Initializes Outcome plot
% GUI plugin displays the settings from the "GUI" subfield of a settings struct.
BpodParameterGUI('init', S); % Initialize parameter GUI plugin--Creates a user interface for viewing and manual override
% Initialize performance graph
performancefigure = figure('Name','OutcomesGraph','NumberTitle','off', 'Position', [1400 30 500 900]); % open appropriate figure
b = categorical({'Hits','Miss','CR', 'FA'}); c = reordercats(b, {'Hits', 'Miss', 'CR', 'FA'}); % set x-axis
subplot(7,1,1);hold on; % graph 1
hit = 0; miss = 0; cr = 0; fa =0; % initalizes numbers of hits, miss, cr, fa to 0
z = [0 0 0 0]; % initiates array of graph
OutcomesGraph = bar(gca,c , z); title('Reinforcement'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:20:110); % Performance figure
numHit = text(1:length(c(1)),z(1),num2str(hit),'HorizontalAlignment','center','VerticalAlignment','bottom'); % initalizes number of hits
numMiss = text(2,z(2), num2str(miss), 'HorizontalAlignment','center','VerticalAlignment','bottom'); % initalizes number of misses
numCR = text(3,z(3), num2str(cr), 'HorizontalAlignment','center','VerticalAlignment','bottom'); % initializes number of correct rejects
numFA = text(4,z(4), num2str(fa), 'HorizontalAlignment','center','VerticalAlignment','bottom'); % initializes number of false alarms
subplot(7,1,2);hold on; % graph 2
hit2 = 0; miss2 = 0; cr2 = 0; fa2 =0; % initializes numbers of hits, miss, cr, fa in probe to 0
z2 = [0 0 0 0]; % initiates array of graph
ProbeGraph = bar(gca, c, z2); title('Probe'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:20:110); % performance in probe
numHit2 = text(1:length(c(1)),z2(1),num2str(hit2),'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 2); % initalizes number of hits
numMiss2 = text(2,z2(2), num2str(miss2), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 2); % initalizes number of misses
numCR2 = text(3, z2(3),num2str(cr2),'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 2); % initializes number of correct rejects
numFA2 = text(4,z2(4), num2str(fa2), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 2); % initializes number of false alarms
subplot(7, 1, 3);hold on; % graph 3
hit3 = 0; miss3 = 0; cr3 = 0; fa3 =0; % initializes numbers of hits, miss, cr, fa in probe to 0
z3 = [0 0 0 0]; % initiates array of graph
FullTrialOptoGraph = bar(gca, c, z3); title('Fifth Dead Opto'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:20:110); % performance in probe
numHit3 = text(1:length(c(1)),z3(1),num2str(hit3),'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 3); % initalizes number of hits
numMiss3 = text(2,z3(2), num2str(miss3), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 3); % initalizes number of misses
numCR3 = text(3, z3(3),num2str(cr3),'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 3); % initializes number of correct rejects
numFA3 = text(4,z3(4), num2str(fa3), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(7, 1, 3); % initializes number of false alarms
na =0;
subplot(7,1,4);hold on; % graph 4
hit4 = 0; miss4 = 0; cr4 = 0; fa4 =0; % initializes numbers of hits, miss, cr, fa in catch light-off to 0
z4 = [0 0 0 0]; % initiates array of graph
CatchOFFGraph = bar(gca, c, z4); title('Catch light-OFF'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:20:110); % performance
numHit4 = text(1:length(c(1)),z4(1),num2str(hit4),'HorizontalAlignment','center','VerticalAlignment','bottom'); % initalizes number of hits
numMiss4 = text(2,z4(2), num2str(miss4), 'HorizontalAlignment','center','VerticalAlignment','bottom'); % initalizes number of misses
numCR4 = text(3, z4(3),num2str(cr4),'HorizontalAlignment','center','VerticalAlignment','bottom'); % initializes number of correct rejects
numFA4 = text(4,z4(4), num2str(fa4), 'HorizontalAlignment','center','VerticalAlignment','bottom'); % initializes number of false alarms
subplot(7,1,5);hold on; % graph 5
hit5 = 0; miss5 = 0; cr5 = 0; fa5 =0; % initializes numbers of hits, miss, cr, fa in catch light-on to 0
z5 = [0 0 0 0]; % initiates array of graph
CatchONGraph = bar(gca, c, z5); title('Catch light-ON'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:20:110); % performance
numHit5 = text(1:length(c(1)),z5(1),num2str(hit5),'HorizontalAlignment','center','VerticalAlignment','bottom'); % initalizes number of hits
numMiss5 = text(2,z5(2), num2str(miss5), 'HorizontalAlignment','center','VerticalAlignment','bottom'); % initalizes number of misses
numCR5 = text(3, z5(3),num2str(cr5),'HorizontalAlignment','center','VerticalAlignment','bottom'); % initializes number of correct rejects
numFA5 = text(4,z5(4), num2str(fa5), 'HorizontalAlignment','center','VerticalAlignment','bottom'); % initializes number of false alarms
na =0;
subplot(7,1,6);hold on; % graph 6 Tone Opto
hit6 = 0; miss6 = 0; cr6 = 0; fa6 =0; % initializes numbers of hits, miss, cr, fa in catch light-on to 0
z6 = [0 0 0 0]; % initiates array of graph
ToneOptoGraph = bar(gca, c, z6); title('Second Dead Opto'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:20:110); % performance
numHit6 = text(1:length(c(1)),z6(1),num2str(hit5),'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 6); % initalizes number of hits
numMiss6 = text(2,z6(2), num2str(miss6), 'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 6); % initalizes number of misses
numCR6 = text(3, z6(3),num2str(cr6),'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 6); % initializes number of correct rejects
numFA6 = text(4,z6(4), num2str(fa6), 'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 6); % initializes number of false alarms
na =0;
subplot(7,1,7);hold on; % graph 7 Chocie Opto
hit7 = 0; miss7 = 0; cr7 = 0; fa7 =0; % initializes numbers of hits, miss, cr, fa in catch light-on to 0
z7 = [0 0 0 0]; % initiates array of graph
ChoiceOptoGraph = bar(gca, c, z7); title('First Dead Opto'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:20:110); % performance
numHit7 = text(1:length(c(1)),z7(1),num2str(hit7),'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 7); % initalizes number of hits
numMiss7 = text(2,z7(2), num2str(miss7), 'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 7); % initalizes number of misses
numCR7 = text(3, z7(3),num2str(cr7),'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 7); % initializes number of correct rejects
numFA7 = text(4,z7(4), num2str(fa7), 'HorizontalAlignment','center','VerticalAlignment','bottom');subplot(7, 1, 7); % initializes number of false alarms
na =0;
%% Define stimuli and send to sound server
S.SF = 192000; % Sound card sampling rate
% Program sound server
% PsychToolboxSoundServer('init')
if ~isfield(BpodSystem.PluginObjects, 'Sound')
    BpodSystem.PluginObjects.Sound = PsychToolboxAudio;
end
GoFreq = go_factor*(GenerateSineWave(S.SF, S.GUI.SinWaveFreqGo, S.GUI.SoundDuration)); % Sampling freq (hz), Sine frequency (hz), duration (s)
NoGoFreq = nogo_factor*(GenerateSineWave(S.SF, S.GUI.SinWaveFreqNoGo, S.GUI.SoundDuration)); % Sampling freq (hz), Sine frequency (hz), duration (s)

BpodSystem.PluginObjects.Sound.load(1, GoFreq); % Load specified sound within trial
BpodSystem.PluginObjects.Sound.load(2, NoGoFreq); % Load specified sound within trial

% Set soft code handler to trigger sounds
BpodSystem.SoftCodeHandlerFunction = 'SoftCodeHandler_PlaySound';
sma = PrepareStateMachine(S, TrialTypes, 1, []); % Prepare state machine for trial 1 with empty "current events" variable
TrialManager.startTrial(sma); % Sends & starts running first trial's state machine. A MATLAB timer object updates the 
                              % console UI, while code below proceeds in parallel.
% In this case, we don't need trial events to build the state machine - but
% they are available in currentTrialEvents.
%% Main trial loop
for currentTrial = 1:ntrials
    currentTrialEvents = TrialManager.getCurrentEvents({'WaitForLick', 'OpenValve'}); % Hangs here until Bpod enters one of the listed trigger states, then returns current trial's states visited + events captured to this point
    if BpodSystem.Status.BeingUsed == 0; return; end % If user hit console "stop" button, end session 
    [sma, S] = PrepareStateMachine(S, TrialTypes, currentTrial+1, currentTrialEvents); % Prepare next state machine.
    % Since PrepareStateMachine is a function with a separate workspace, pass any local variables needed to make 
    % the state machine as fields of settings struct S e.g. S.learningRate = 0.2.
    SendStateMachine(sma, 'RunASAP'); % With TrialManager, you can send the next trial's state machine while the current trial is ongoing
    RawEvents = TrialManager.getTrialData; % Hangs here until trial is over, then retrieves full trial's raw data
    if BpodSystem.Status.BeingUsed == 0; return; end % If user hit console "stop" button, end session 
    HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
    TrialManager.startTrial(sma); % Start processing the next trial's events (** can call with no argument since SM was already sent)
    if ~isempty(fieldnames(RawEvents)) % If trial data was returned from last trial, update plots and save data
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
        BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
        BpodSystem.Data.TrialTypes(currentTrial) = TrialTypes(currentTrial); % Adds the trial type of the current trial to data
        UpdateOutcomePlot(TrialTypes, BpodSystem.Data);
        SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
    end  
    % Updates performance graphs
    set(performancefigure,'Name',['Next trial: ',num2str(currentTrial+1)]);
    if S.context(currentTrial) == 0 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.OpenValve) % hit
        hit2 = hit2+1; % updates number of hits in probe
        z2=[((hit2/(hit2+miss2))*100) ((miss2/(hit2+miss2))*100) z2(3) z2(4)]; % change percentage difference between hit and miss
        set(ProbeGraph, 'YData', z2); % Updates probe performance plot
    elseif S.context(currentTrial) == 0 & TrialTypes(currentTrial) == 3 % go
        miss2 = miss2 +1; % update num of misses in probe
        z2=[((hit2/(hit2+miss2))*100) ((miss2/(hit2+miss2))*100) z2(3) z2(4)]; % change percentage difference between hit and miss
        set(ProbeGraph, 'YData', z2); % Update probe performance plot
    elseif S.context(currentTrial) == 0 & TrialTypes(currentTrial) == 4 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject) % If no-go and correct reject
        cr2 = cr2 +1; % update num of correct rejects
        z2 = [z2(1) z2(2) ((cr2/(cr2+fa2))*100) ((fa2/(fa2+cr2))*100)]; % change percentage difference between false alarm and correct rejects
        set(ProbeGraph, 'YData', z2); % Update performance plot
    elseif S.context(currentTrial) == 0 & TrialTypes(currentTrial) == 4 % no go
        fa2 = fa2+1; % update num of false alarms
        z2 = [z2(1) z2(2) ((cr2/(cr2+fa2))*100) ((fa2/(fa2+cr2))*100)]; % change percentage difference between false alarm and correct rejects
        set(ProbeGraph, 'YData', z2); % Update performance plot
        
        % reinforced
    elseif S.context(currentTrial) == 2 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.OpenValve)
        hit = hit +1; % update num of hits
        z=[((hit/(hit+miss))*100) ((miss/(hit+miss))*100) z(3) z(4)]; % change percentage difference between hit and miss
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
    elseif S.context(currentTrial) == 2 & TrialTypes(currentTrial) == 1 % go
        miss = miss +1; % update num of misses
        z=[((hit/(hit+miss))*100) ((miss/(hit+miss))*100) z(3) z(4)]; % change percentage difference between hit and miss
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
    elseif S.context(currentTrial) == 2 & TrialTypes(currentTrial) == 2 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject) % If no-go and correct reject
        cr = cr +1; % update num of correct rejects
        z = [z(1) z(2) ((cr/(cr+fa))*100) ((fa/(fa+cr))*100)]; % change percentage difference between false alarm and correct rejects
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
    elseif S.context(currentTrial) == 2 & TrialTypes(currentTrial) == 2 % no go
        fa = fa+1; % update num of false alarms
        z = [z(1) z(2) ((cr/(cr+fa))*100) ((fa/(fa+cr))*100)]; % change percentage difference between false alarm and correct rejects
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
        
        % full trial opto reinforced
    elseif S.context(currentTrial) == 1 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.OpenValve)
        hit3 = hit3 +1; % update num of hits
        z3=[((hit3/(hit3+miss3))*100) ((miss3/(hit3+miss3))*100) z3(3) z3(4)]; % change percentage difference between hit and miss
        set(FullTrialOptoGraph, 'YData', z3); % Update opto performance plot
    elseif S.context(currentTrial) == 1 & TrialTypes(currentTrial) == 5 % go
        miss3 = miss3 +1; % update num of misses
        z3=[((hit3/(hit3+miss3))*100) ((miss3/(hit3+miss3))*100) z3(3) z3(4)]; % change percentage difference between hit and miss
        set(FullTrialOptoGraph, 'YData', z3); % Update opto performance plot
    elseif S.context(currentTrial) == 1 & TrialTypes(currentTrial) == 6 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject) % If no-go and correct reject
        cr3 = cr3 +1; % update num of correct rejects
        z3 = [z3(1) z3(2) ((cr3/(cr3+fa3))*100) ((fa3/(fa3+cr3))*100)]; % change percentage difference between false alarm and correct rejects
        set(FullTrialOptoGraph, 'YData', z3); % Update opto performance plot
    elseif S.context(currentTrial) == 1 & TrialTypes(currentTrial) == 6 % no go
        fa3 = fa3+1; % update num of false alarms
        z3 = [z3(1) z3(2) ((cr3/(cr3+fa3))*100) ((fa3/(fa3+cr3))*100)]; % change percentage difference between false alarm and correct rejects
        set(FullTrialOptoGraph, 'YData', z3); % Update opto performance plot
   
    elseif (S.context(currentTrial) == 3 & TrialTypes(currentTrial) == 7 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject)) ... % no go
            | (S.context(currentTrial) == 3 & TrialTypes(currentTrial) == 8 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject))
        cr4 = cr4 +1; % update num of correct rejects
        z4 = [z4(1) z4(2) ((cr4/(cr4+fa4))*100) ((fa4/(fa4+cr4))*100)]; % change percentage difference between false alarm and correct rejects
        set(CatchOFFGraph, 'YData', z4); % Update opto performance plot
    elseif (S.context(currentTrial) == 4 & TrialTypes(currentTrial) == 9 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject)) ... % no go
            | (S.context(currentTrial) == 4 & TrialTypes(currentTrial) == 10 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject))
        cr5 = cr5+1; % update num of false alarms
        z5 = [z5(1) z5(2) ((cr5/(cr5+fa5))*100) ((fa5/(fa5+cr5))*100)]; % change percentage difference between false alarm and correct rejects
        set(CatchONGraph, 'YData', z5); % Update opto performance plot
    elseif (S.context(currentTrial) == 3 & TrialTypes(currentTrial) == 7) ... % no go
            | (S.context(currentTrial) == 3 & TrialTypes(currentTrial) == 8)
        fa4 = fa4+1; % update num of false alarms
        z4 = [z4(1) z4(2) ((cr4/(cr4+fa4))*100) ((fa4/(fa4+cr4))*100)]; % change percentage difference between false alarm and correct rejects
        set(CatchOFFGraph, 'YData', z4); % Update opto performance plot
    elseif (S.context(currentTrial) == 4 & TrialTypes(currentTrial) == 9) ... % no go
            | (S.context(currentTrial) == 4 & TrialTypes(currentTrial) == 10)
        fa5 = fa5+1; % update num of false alarms
        z5 = [z5(1) z5(2) ((cr5/(cr5+fa5))*100) ((fa5/(fa5+cr5))*100)]; % change percentage difference between false alarm and correct rejects
        set(CatchONGraph, 'YData', z5); % Update opto performance plot        
    elseif S.context(currentTrial) == 3
        na = na+1;
        
    elseif S.context(currentTrial) == 5 & TrialTypes(currentTrial) == 11 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.OpenValve)% Tone Opto go
        hit6 = hit6 +1; % update num of hits
        z6=[((hit6/(hit6+miss6))*100) ((miss6/(hit6+miss6))*100) z6(3) z6(4)]; % change percentage difference between hit and miss
        set(ToneOptoGraph, 'YData', z6); % Update opto performance plot
    elseif S.context(currentTrial) == 5 & TrialTypes(currentTrial) == 11 % Tone Opto go
        miss6 = miss6 +1; % update num of misses
        z6=[((hit6/(hit6+miss6))*100) ((miss6/(hit6+miss6))*100) z6(3) z6(4)]; % change percentage difference between hit and miss
        set(ToneOptoGraph, 'YData', z6); % Update opto performance plot
    elseif S.context(currentTrial) == 5 & TrialTypes(currentTrial) == 12 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject) % If no-go and correct reject
        cr6 = cr6 +1; % update num of correct rejects
        z6 = [z6(1) z6(2) ((cr6/(cr6+fa6))*100) ((fa6/(fa6+cr6))*100)]; % change percentage difference between false alarm and correct rejects
        set(ToneOptoGraph, 'YData', z6); % Update opto performance plot
    elseif (S.context(currentTrial) == 5 & TrialTypes(currentTrial) == 12) ... % Tone Opto no go
        fa6 = fa6+1; % update num of false alarms
        z6 = [z6(1) z6(2) ((cr6/(cr6+fa6))*100) ((fa6/(fa6+cr6))*100)]; % change percentage difference between false alarm and correct rejects
        set(ToneOptoGraph, 'YData', z6); % Update Tone opto performance plot 
        
        
     elseif S.context(currentTrial) == 6 & TrialTypes(currentTrial) == 13 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.OpenValve)% Choice Opto go
        hit7 = hit7 +1; % update num of hits
        z7=[((hit7/(hit7+miss7))*100) ((miss7/(hit7+miss7))*100) z7(3) z7(4)]; % change percentage difference between hit and miss
        set(ChoiceOptoGraph, 'YData', z7); % Update opto performance plot
    elseif S.context(currentTrial) == 6 & TrialTypes(currentTrial) == 13 % Choice Opto go
        miss7 = miss7 +1; % update num of misses
        z7=[((hit7/(hit7+miss7))*100) ((miss7/(hit7+miss7))*100) z7(3) z7(4)]; % change percentage difference between hit and miss
        set(ChoiceOptoGraph, 'YData', z7); % Update opto performance plot
    elseif S.context(currentTrial) == 6 & TrialTypes(currentTrial) == 14 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject) % If no-go and correct reject
        cr7 = cr7 +1; % update num of correct rejects
        z7 = [z7(1) z7(2) ((cr7/(cr7+fa7))*100) ((fa7/(fa7+cr7))*100)]; % change percentage difference between false alarm and correct rejects
        set(ChoiceOptoGraph, 'YData', z7); % Update opto performance plot
    elseif (S.context(currentTrial) == 6 & TrialTypes(currentTrial) == 14) ... % Tone Opto no go
        fa7 = fa7+1; % update num of false alarms
        z7 = [z7(1) z7(2) ((cr7/(cr7+fa7))*100) ((fa7/(fa7+cr7))*100)]; % change percentage difference between false alarm and correct rejects
        set(ChoiceOptoGraph, 'YData', z7); % Update Tone opto performance plot 
    end
% Update numbers of hit, misses, cr, and false alarms above each bin in performance graphs (deletes previous, updates number in specific plot)    
delete(numHit2); subplot(712); numHit2 = text(1:length(c(1)),z2(1),num2str(hit2),'vert','bottom','horiz','center'); 
delete(numMiss2);subplot(712); numMiss2 = text(2,z2(2), num2str(miss2), 'vert', 'bottom', 'horiz', 'center');
delete(numCR2);subplot(712);numCR2 = text(3,z2(3), num2str(cr2), 'vert', 'bottom', 'horiz', 'center');
delete(numFA2);subplot(712);numFA2 = text(4,z2(4), num2str(fa2), 'vert', 'bottom', 'horiz', 'center');
delete(numHit);subplot(711);numHit = text(1:length(c(1)),z(1),num2str(hit),'HorizontalAlignment','center','VerticalAlignment','bottom');
delete(numMiss);subplot(711);numMiss = text(2,z(2), num2str(miss), 'vert', 'bottom', 'horiz', 'center');
delete(numCR);subplot(711);numCR = text(3,z(3), num2str(cr), 'vert', 'bottom', 'horiz', 'center');
delete(numFA);subplot(711);numFA = text(4,z(4), num2str(fa), 'vert', 'bottom', 'horiz', 'center');
delete(numHit3);subplot(713);numHit3 = text(1:length(c(1)),z3(1),num2str(hit3),'HorizontalAlignment','center','VerticalAlignment','bottom');
delete(numMiss3);subplot(713);numMiss3 = text(2,z3(2), num2str(miss3), 'vert', 'bottom', 'horiz', 'center');
delete(numCR3);subplot(713);numCR3 = text(3,z3(3), num2str(cr3), 'vert', 'bottom', 'horiz', 'center');
delete(numFA3);subplot(713);numFA3 = text(4,z3(4), num2str(fa3), 'vert', 'bottom', 'horiz', 'center');
delete(numHit4);subplot(714);numHit4 = text(1:length(c(1)),z4(1),num2str(hit4),'HorizontalAlignment','center','VerticalAlignment','bottom');
delete(numMiss4);subplot(714);numMiss4 = text(2,z4(2), num2str(miss4), 'vert', 'bottom', 'horiz', 'center');
delete(numCR4);subplot(714);numCR4 = text(3,z4(3), num2str(cr4), 'vert', 'bottom', 'horiz', 'center');
delete(numFA4);subplot(714);numFA4 = text(4,z4(4), num2str(fa4), 'vert', 'bottom', 'horiz', 'center');
delete(numHit5);subplot(715);numHit5 = text(1:length(c(1)),z5(1),num2str(hit5),'HorizontalAlignment','center','VerticalAlignment','bottom');
delete(numMiss5);subplot(715);numMiss5 = text(2,z5(2), num2str(miss5), 'vert', 'bottom', 'horiz', 'center');
delete(numCR5);subplot(715);numCR5 = text(3,z5(3), num2str(cr5), 'vert', 'bottom', 'horiz', 'center');
delete(numFA5);subplot(715);numFA5 = text(4,z5(4), num2str(fa5), 'vert', 'bottom', 'horiz', 'center');

delete(numHit6);subplot(716);numHit6 = text(1:length(c(1)),z6(1),num2str(hit6),'HorizontalAlignment','center','VerticalAlignment','bottom');
delete(numMiss6);subplot(716);numMiss6 = text(2,z6(2), num2str(miss6), 'vert', 'bottom', 'horiz', 'center');
delete(numCR6);subplot(716);numCR6 = text(3,z6(3), num2str(cr6), 'vert', 'bottom', 'horiz', 'center');
delete(numFA6);subplot(716);numFA6 = text(4,z6(4), num2str(fa6), 'vert', 'bottom', 'horiz', 'center');
delete(numHit7);subplot(717);numHit7 = text(1:length(c(1)),z7(1),num2str(hit7),'HorizontalAlignment','center','VerticalAlignment','bottom');
delete(numMiss7);subplot(717);numMiss7 = text(2,z7(2), num2str(miss7), 'vert', 'bottom', 'horiz', 'center');
delete(numCR7);subplot(717);numCR7 = text(3,z7(3), num2str(cr7), 'vert', 'bottom', 'horiz', 'center');
delete(numFA7);subplot(717);numFA7 = text(4,z7(4), num2str(fa7), 'vert', 'bottom', 'horiz', 'center');
end
function [sma, S] = PrepareStateMachine(S, TrialTypes, currentTrial, ~)

% Default values
hitdelay = 4;
missdelay = 2;
crdelay = 2;
fadelay = 7;
deadperiod = 0.2;
smalldeadperiod = 0.1;
smallestdeadperiod = 0.05;
smalldeadperiod2=0.075;
waitfornolick = 1;
respwindow = 2.5;
hitdelay1 = 2.5;
hitdelay2 = 1.5;
fadelay1 = 2.5;
fadelay2 = 4.5;
lightOnBeforeTone = 0.1;
lightOnDuringHit = 2;
lightOffDuringHit = 2;
if lightOnDuringHit+lightOffDuringHit ~= hitdelay
    error(['Hit delay should always be equal to: ' num2str(hitdelay)]);
end
lightOnDuringFadelay = 2;
lightOffDuringFadelay = 5;
if lightOnDuringFadelay+lightOffDuringFadelay ~= fadelay
    error(['FA delay should always be equal to: ' num2str(fadelay)]);
end
delaypostlick = 0.07;
lightonpostlick = 0.1;
sma = NewStateMatrix(); % Assemble state matrix
S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin
R = GetValveTimes(S.GUI.RewardAmount, 1); ValveTime = R;  % Update reward amounts
switch TrialTypes(currentTrial) % Determine trial-specific state matrix fields
    case 1 % GO trial
        Stimulus = 1; % GoTone
        INResponse = 'OpenValve';
        NOResponse = 'Miss';
    case 2 % No-Go trial
        Stimulus = 2; %NoGoTone
        INResponse = 'Punish';
        NOResponse = 'CorrectReject'; 
    case 3 % Probe-GO trial
        Stimulus = 1; % GoTone
        INResponse = 'OpenValve';
        NOResponse = 'Miss';
    case 4 % Probe-No-Go trial  
        Stimulus = 2; % NoGoTone
        INResponse = 'Punish';
        NOResponse = 'CorrectReject'; 
    case 5 % Dead Third Opto-GO trial
        Stimulus = 1; % GoTone
        INResponse = 'OpenValve';
        NOResponse = 'Miss';
        LightON = 1;
    case 6 % Dead Third Opto-No-Go trial
        Stimulus = 2; %NoGoTone
        INResponse = 'Punish';
        NOResponse = 'CorrectReject'; 
        LightON = 1;
    case 7 % Catch reinforced light OFF (No target)
        INResponse = 'Punish';
        NOResponse = 'CorrectReject'; 
    case 8 % Catch reinforced light OFF (No foil)
        INResponse = 'Punish';
        NOResponse = 'CorrectReject'; 
    case 9 % Catch reinforced light ON (No target)
        INResponse = 'PunishLightON';
        NOResponse = 'CorrectReject'; 
        LightON = 1;
    case 10 % Catch reinforced light ON (No foil)
        INResponse = 'PunishLightON';
        NOResponse = 'CorrectReject'; 
        LightON = 1;    
    case 11 % Second Dead Opto-GO trial
        Stimulus = 1; % GoTone
        INResponse = 'OpenValve';
        NOResponse = 'Miss';
        LightON = 1;
    case 12 % Second Dead Opto-No-Go trial
        Stimulus = 2; %NoGoTone
        INResponse = 'Punish';
        NOResponse = 'CorrectReject'; 
        LightON = 1;
    case 13 % First Dead Opto-GO trial
        Stimulus = 1; % GoTone
        INResponse = 'OpenValve';
        NOResponse = 'Miss';
        LightON = 1;
    case 14 % First Dead Opto-No-Go trial
        Stimulus = 2; %NoGoTone
        INResponse = 'PunishFirst';
        NOResponse = 'CorrectReject'; 
        LightON = 1;
end
if S.context(currentTrial) == 2 % Reinforced context 
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', waitfornolick, ...
        'StateChangeConditions', {'Tup', 'Stimulus', 'Port2In', 'Stop'}, ... % if no action, stimulus activated; if action, stop period
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stop', ... % Stop period to assure no activity during Pre-Trial
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'Dead'}, ...
        'OutputActions', {'SoftCode', Stimulus});  
    sma = AddState(sma, 'Name', 'Dead', ... % 100ms dead period
        'Timer', deadperiod, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', respwindow, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 1}); 
    sma = AddState(sma, 'Name', 'Drinking', ... % 4 seconds for drinking
        'Timer', hitdelay,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'ValveState', 2}); 
    sma = AddState(sma, 'Name', 'Miss', ... % 2 second time out for miss
        'Timer', missdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', crdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Punish', ... % 7 second punish state
        'Timer', fadelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {}); 
    
elseif S.context(currentTrial) == 1 % Fifth Dead context 
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', waitfornolick, ...
        'StateChangeConditions', {'Tup', 'Stimulus','Port2In', 'Stop'}, ... % if no action, trigger light ON; if action, stop period
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stop', ... 
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone presentation
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'SmallFirstDead'}, ...
        'OutputActions', {'SoftCode', Stimulus});  
    sma = AddState(sma, 'Name', 'SmallFirstDead', ... % 75ms dead period
        'Timer', smalldeadperiod2, ...
        'StateChangeConditions', {'Tup', 'MiddleDead'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'MiddleDead', ... % 50ms dead period
        'Timer', smallestdeadperiod, ...
        'StateChangeConditions', {'Tup', 'LastDead'}, ...
        'OutputActions', {'BNC1', LightON}); 
    sma = AddState(sma, 'Name', 'LastDead', ... % 75ms dead period
        'Timer', smalldeadperiod2, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {});     
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', respwindow, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {});    
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 1, }); 
    sma = AddState(sma, 'Name', 'Drinking', ... % 4 seconds for drinking
        'Timer', hitdelay,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'ValveState', 2}); 
    sma = AddState(sma, 'Name', 'Miss', ... % 2 second time out for miss
        'Timer', missdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', crdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Punish', ... % 7 second punish state
        'Timer', fadelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {}); 
    
elseif S.context(currentTrial) == 3 % CATCH context, reinforced NO light 
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', waitfornolick, ...
        'StateChangeConditions', {'Tup', 'Stimulus', 'Port2In', 'Stop'}, ... % if no action, stimulus activated; if action, stop period
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stop', ... % Stop period to assure no activity during Pre-Trial
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stimulus', ... % No tone
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'Dead'}, ...
        'OutputActions', {});  
    sma = AddState(sma, 'Name', 'Dead', ... % 100ms dead period
        'Timer', deadperiod, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', respwindow, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, FA. If no lick, CR period
        'OutputActions', {});     
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', crdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Punish', ... % 7 second punish state
        'Timer', fadelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {}); 
    
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 1}); 
    sma = AddState(sma, 'Name', 'Drinking', ... % 4 seconds for drinking
        'Timer', hitdelay,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'ValveState', 2}); 
    
elseif S.context(currentTrial) == 4 % CATCH context, reinforced light ON
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', waitfornolick, ...
        'StateChangeConditions', {'Tup', 'LightON', 'Port2In', 'Stop'}, ... % if no action, trigger light ON; if action, stop period
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stop', ... 
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'LightON', ... % Turn the light ON during 0.1s before tone presentation
        'Timer', lightOnBeforeTone, ...
        'StateChangeConditions', {'Tup', 'Stimulus'}, ... 
        'OutputActions', {'BNC1', LightON}); 
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone presentation
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'Dead'}, ...
        'OutputActions', {'BNC1', LightON});  
    sma = AddState(sma, 'Name', 'Dead', ... % 200ms dead period
        'Timer', deadperiod, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {'BNC1', LightON}); 
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', respwindow, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {'BNC1', LightON});    
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', crdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {});     
    sma = AddState(sma, 'Name', 'PunishLightON', ... % % 1st second of punish state with light
        'Timer', lightOnDuringFadelay,...
        'StateChangeConditions', {'Tup', 'PunishLightOFF'}, ...
        'OutputActions', {'BNC1', LightON}); 
    sma = AddState(sma, 'Name', 'PunishLightOFF', ... % 6 more seconds of punish state without light
        'Timer', lightOffDuringFadelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {});     
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {}); 
    
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 1}); 
    sma = AddState(sma, 'Name', 'Drinking', ... % 4 seconds for drinking
        'Timer', hitdelay,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'ValveState', 2}); 
    
% context 5, light on only during second dead period (100ms)
elseif S.context(currentTrial) == 5 % second dead OPTO context 
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', waitfornolick, ...
        'StateChangeConditions', {'Tup', 'LightON', 'Port2In', 'Stop'}, ... % if no action, trigger light ON; if action, stop period
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stop', ... 
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'LightON', ...
        'Timer', lightOnBeforeTone, ...
        'StateChangeConditions', {'Tup', 'Stimulus'}, ... 
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone presentation
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'DeadFirst'}, ...
        'OutputActions', {'SoftCode', Stimulus});  
    sma = AddState(sma, 'Name', 'DeadFirst', ... % 100ms dead period 
        'Timer', smalldeadperiod, ...
        'StateChangeConditions', {'Tup', 'DeadSecond'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'DeadSecond', ... % 100ms dead period turn the light on
        'Timer', smalldeadperiod, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {'BNC1', LightON}); 
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', respwindow, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {'BNC1', LightON});    
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'DrinkingLightOFF'},...
        'OutputActions', {'ValveState', 1}); 
    sma = AddState(sma, 'Name', 'DrinkingLightON', ... % 1st second of drinking with light
        'Timer', lightOnDuringHit,...
        'StateChangeConditions', {'Tup', 'DrinkingLightOFF'},...
        'OutputActions', {'ValveState', 2}); 
    sma = AddState(sma, 'Name', 'DrinkingLightOFF', ... % 3 more seconds for drinking without light
        'Timer', lightOffDuringHit,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'ValveState', 2}); 
    
    sma = AddState(sma, 'Name', 'Miss', ... % 2 second time out for miss
        'Timer', missdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', crdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Punish', ... % 7 second punish state
        'Timer', fadelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {});    
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {}); 
    
% context 6, light on during first 100ms period in delay
elseif S.context(currentTrial) == 6 % first 100ms delay OPTO context 
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', waitfornolick, ...
        'StateChangeConditions', {'Tup', 'Stimulus', 'Port2In', 'Stop'}, ... % if no action, stimulus activated; if action, stop period
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stop', ... % Stop period to assure no activity during Pre-Trial
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'DeadOn'}, ...
        'OutputActions', {'SoftCode', Stimulus});  
    sma = AddState(sma, 'Name', 'DeadOn', ... % 100ms dead period
        'Timer', smalldeadperiod, ...
        'StateChangeConditions', {'Tup', 'DeadOff'}, ...
        'OutputActions', {'BNC1', LightON});
    sma = AddState(sma, 'Name', 'DeadOff', ... % 100ms dead period, no light
        'Timer', smalldeadperiod, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', respwindow, ... 
        'StateChangeConditions', {'Port1In','DelayPostOperantLick', 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {});  
    sma = AddState(sma, 'Name', 'DelayPostOperantLick', ... % Open valve for reward
        'Timer', delaypostlick,... % this only exists in celine's post lick code? 2/9/2024
        'StateChangeConditions', {'Tup', INResponse},...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'DrinkingFirst'},...
        'OutputActions', {'ValveState', 1}); 
    sma = AddState(sma, 'Name', 'DrinkingFirst', ... % 2.5 seconds for drinking with light
        'Timer', hitdelay1,...
        'StateChangeConditions', {'Tup', 'DrinkingSecond'},...
        'OutputActions', {'ValveState', 2}); 
    sma = AddState(sma, 'Name', 'DrinkingSecond', ... % 1.5 seconds for drinking without light
        'Timer', hitdelay2,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'ValveState', 2}); 
    sma = AddState(sma, 'Name', 'Miss', ... % 2 second time out for miss
        'Timer', missdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', crdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'PunishFirst', ... % 7 second punish state
        'Timer', fadelay1,...
        'StateChangeConditions', {'Tup', 'PunishSecond'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'PunishSecond', ... % 7 second punish state
        'Timer', fadelay2,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {}); 
    
elseif S.context(currentTrial) == 0 % Probe context
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', waitfornolick, ...
        'StateChangeConditions', {'Tup', 'Stimulus', 'Port1In', 'Stop'}, ... % if no action, stimulus activated; if action, stop period
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Stop', ... % Stop period to assure no activity during Pre-Trial
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'Dead'}, ...
        'OutputActions', {'SoftCode', Stimulus}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'Dead', ... % 100ms dead period
        'Timer', deadperiod, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', respwindow, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Drinking', ... % 4 seconds for drinking
        'Timer', hitdelay,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Miss', ... % 2 second time out for miss
        'Timer', missdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', crdelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Punish', ... % 7 second punish state
        'Timer', fadelay,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {}); % Keeps lick tube out
end

function UpdateOutcomePlot(TrialTypes, Data)
global BpodSystem
Outcomes = zeros(1,Data.nTrials);% Creates a vector for each completed trial, listing outcomes
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.OpenValve)
        Outcomes(x) = -1; % green circle for hits
    elseif TrialTypes(x) == 1 % go
        Outcomes(x) = 1; % green x for Miss
    elseif TrialTypes(x) == 3 % go
        Outcomes(x) = 1; % green x for Miss
    elseif TrialTypes(x) == 5 % go
        Outcomes(x) = 1; % green x for Miss    
    elseif TrialTypes(x) == 2 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject) % If no-go and correct reject
        Outcomes(x) = 2; % unfilled green circle for Correct Rejects
    elseif TrialTypes(x) == 6 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject) % If no-go and correct reject
        Outcomes(x) = 2; % unfilled green circle for Correct Rejects
    elseif TrialTypes(x) == 4 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject) % If no-go and correct reject
        Outcomes(x) = 2; % unfilled green circle for Correct Rejects
    elseif TrialTypes(x) == 2 % no go
        Outcomes(x) = 0; % red X for false alarm
    elseif (TrialTypes(x) == 7 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject))... % If no-go and correct reject 
        | (TrialTypes(x) == 8 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject))...
        | (TrialTypes(x) == 9 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject))...
        | (TrialTypes(x) == 10 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject)) ...
        | (TrialTypes(x) == 14 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject)) 
        Outcomes(x) = 2; % unfilled green circle for Correct Rejects
    elseif TrialTypes(x) == 8 ... % no go
        | TrialTypes(x) == 7 ...
        | TrialTypes(x) == 9 ...
        | TrialTypes(x) == 10 ...
        | TrialTypes(x) == 14
        Outcomes(x) = 3; % red X for false alarm
    end
end
TrialTypeOutcomePlot_Opto(BpodSystem.GUIHandles.OutcomePlot,'update',Data.nTrials+1,TrialTypes',Outcomes);
