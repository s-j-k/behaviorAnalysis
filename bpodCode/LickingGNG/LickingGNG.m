function LickingGNG
% this code makes the lick cup stay down between reward deliveries. 
% updated 20231015 
global BpodSystem
global Amount
%% WHEN NOT USING PROBE, COMMENT OUT LINES 21, 24, 32, 33
% input weight
weight = 28.1;

%% Create trial manager object
TrialManager = TrialManagerObject;

%% Define parameters
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    S.GUI.RewardAmount = 3; % ul % due to calibration for output
    S.GUI.SoundDuration = 0.1; % duration of sound
    S.GUI.SinWaveFreqGo = 11314; % Frequency of go cue
    S.GUI.SinWaveFreqNoGo = 19027; % Frequency of no-go cue
    S.GUI.Weight = weight;
    S.GUIPanels.Sound = {'SinWaveFreqGo', 'SinWaveFreqNoGo', 'SoundDuration'}; % Labels for sound panel
end
Amount = S.GUI.RewardAmount;%actual calibration amount (ul from g)
%% Define trials
MaxTrials = 1000; % max trials
n = 0; % first n trials are GO (Type 1) 
probe1= [821 840]; % Input trials for first probe block %% COMMENT WHEN NOT USING PROBE
S.context = ones(MaxTrials, 1); %1 = reinforced context, licktube in
% despite having separate trial types still require different contexts on top of that because separate state machine needs to run with licktube out (multiple valves cannot be open at once)
S.context(probe1(1):probe1(2)) = 0;% 0 = probe context, licktube out %% COMMENT WHEN NOT USING PROBE
randomize = RandStream('mlfg6331_64'); % creates random stream of numbers for data samplign later
TrialTypes = []; % initalizes empty array for Trialtypes
for i = 1:50 % 50 groups of 20 trials, each 20 trials is balanced
    TrialTypes(i,:) = datasample(randomize, [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2],20,'Replace',false);
end
TrialTypes = TrialTypes'; % transposes array for proper formatting
TrialTypes(1:n) = 1; % overwrites first n trials with specified trial
TrialTypes(probe1(1):probe1(1)+1) = 3; % first 2 trials of probe are GO %% COMMENT WHEN NOT USING PROBE
TrialTypes((probe1(1)+2):probe1(2)) = datasample(randomize, [3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4],18,'Replace',false); % balances remaining 18 trials of probe %% COMMENT WHEN NOT USING PROBE
% TrialTypes((probe1(1)+2):probe1(2)) = datasample(randomize, [3 3 3 4 4 4 4 4],8,'Replace',false); % FOR 10 TRIAL PROBE. balances remaining 18 trials of probe %% COMMENT WHEN NOT USING PROBE
BpodSystem.Data.TrialTypes = []; % The trial type of each trial completed will be added here.
tic; % starts timer
%% Initialize plots
BpodSystem.ProtocolFigures.OutcomePlotFig = figure('Position', [10 750 1900 400],'name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off'); % Initializes figure for Outcome plot
BpodSystem.GUIHandles.OutcomePlot = axes('Position', [.075 .3 .9 .6]); % Initializes axes for Outcome plot
% BpodSystem.GUIHandles.ProbeContextLine2 = line((probe2),[-3,-3],'LineStyle','-', 'LineWidth', 10,'Color','b', 'MarkerSize',100); % draws line for probe trials
TrialTypeOutcomePlot(BpodSystem.GUIHandles.OutcomePlot,'init',TrialTypes) % 'ntrials',MaxTrials); % Initializes Outcome plot
% GUI plugin displays the settings from the "GUI" subfield of a settings struct.
BpodParameterGUI('init', S); % Initialize parameter GUI plugin--Creates a user interface for viewing and manual override
% Initialize performance graph
figure('Name','OutcomesGraph','NumberTitle','off', 'Position', [1400 30 500 700]); % open appropriate figure
b = categorical({'Hits','Miss','CR', 'FA'}); c = reordercats(b, {'Hits', 'Miss', 'CR', 'FA'}); % set x-axis
subplot(2,1,1); % graph 1
hit = 0; miss = 0; cr = 0; fa =0; % initalizes numbers of hits, miss, cr, fa to 0
z = [0 0 0 0]; % initiates array of graph
OutcomesGraph = bar(gca,c , z); title('Reinforcement'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:10:110); % Performance figure
numHit = text(1:length(c(1)),z(1),num2str(hit),'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2,1,1); % initalizes number of hits
numMiss = text(2,z(2), num2str(miss), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2,1,1); % initalizes number of misses
numCR = text(3,z(3), num2str(cr), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2,1,1); % initializes number of correct rejects
numFA = text(4,z(4), num2str(fa), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2,1,1);% initializes number of false alarms
subplot(2, 1, 2); % graph 2
hit2 = 0; miss2 = 0; cr2 = 0; fa2 =0; % initializes numbers of hits, miss, cr, fa in probe to 0
z2 = [0 0 0 0]; % initiates array of graph
ProbeGraph = bar(gca, c, z2); title('Probe'); xlabel('Outcome'); ylabel('% Correct'); ylim([0 110]); yticks(0:10:110); % performance in probe
numHit2 = text(1:length(c(1)),z2(1),num2str(hit2),'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2, 1, 2); % initalizes number of hits
numMiss2 = text(2,z2(2), num2str(miss2), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2, 1, 2); % initalizes number of misses
numCR2 = text(3, z2(3),num2str(cr2),'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2, 1, 2); % initializes number of correct rejects
numFA2 = text(4,z2(4), num2str(fa2), 'HorizontalAlignment','center','VerticalAlignment','bottom'); subplot(2, 1, 2); % initializes number of false alarms
%% Define stimuli and send to sound server
S.SF = 192000; % Sound card sampling rate
% Program sound server
if ~isfield(BpodSystem.PluginObjects, 'Sound')
    BpodSystem.PluginObjects.Sound = PsychToolboxAudio;
end
GoFreq = .4*(GenerateSineWave(S.SF, S.GUI.SinWaveFreqGo, S.GUI.SoundDuration)); % Sampling freq (hz), Sine frequency (hz), duration (s)
NoGoFreq = .3*(GenerateSineWave(S.SF, S.GUI.SinWaveFreqNoGo, S.GUI.SoundDuration)); % Sampling freq (hz), Sine frequency (hz), duration (s)

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
for currentTrial = 1:MaxTrials
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
    elseif S.context(currentTrial) == 1 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.OpenValve)
        hit = hit +1; % update num of hits
        z=[((hit/(hit+miss))*100) ((miss/(hit+miss))*100) z(3) z(4)]; % change percentage difference between hit and miss
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
    elseif S.context(currentTrial) == 1 & TrialTypes(currentTrial) == 1 % go
        miss = miss +1; % update num of misses
        z=[((hit/(hit+miss))*100) ((miss/(hit+miss))*100) z(3) z(4)]; % change percentage difference between hit and miss
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
    elseif S.context(currentTrial) == 1 & TrialTypes(currentTrial) == 2 & ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.CorrectReject) % If no-go and correct reject
        cr = cr +1; % update num of correct rejects
        z = [z(1) z(2) ((cr/(cr+fa))*100) ((fa/(fa+cr))*100)]; % change percentage difference between false alarm and correct rejects
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
    elseif S.context(currentTrial) == 1 & TrialTypes(currentTrial) == 2 % no go
        fa = fa+1; % update num of false alarms
        z = [z(1) z(2) ((cr/(cr+fa))*100) ((fa/(fa+cr))*100)]; % change percentage difference between false alarm and correct rejects
        set(OutcomesGraph, 'YData', z); % Update reinforcement performance plot
    end
% Update numbers of hit, misses, cr, and false alarms above each bin in performance graphs (deletes previous, updates number in specific plot)    
delete(numHit2); subplot(212); numHit2 = text(1:length(c(1)),z2(1),num2str(hit2),'vert','bottom','horiz','center'); 
delete(numMiss2);subplot(212); numMiss2 = text(2,z2(2), num2str(miss2), 'vert', 'bottom', 'horiz', 'center');
delete(numCR2);subplot(212);numCR2 = text(3,z2(3), num2str(cr2), 'vert', 'bottom', 'horiz', 'center');
delete(numFA2);subplot(212);numFA2 = text(4,z2(4), num2str(fa2), 'vert', 'bottom', 'horiz', 'center');
delete(numHit);subplot(211);numHit = text(1:length(c(1)),z(1),num2str(hit),'HorizontalAlignment','center','VerticalAlignment','bottom');
delete(numMiss);subplot(211);numMiss = text(2,z(2), num2str(miss), 'vert', 'bottom', 'horiz', 'center');
delete(numCR);subplot(211);numCR = text(3,z(3), num2str(cr), 'vert', 'bottom', 'horiz', 'center');
delete(numFA);subplot(211);numFA = text(4,z(4), num2str(fa), 'vert', 'bottom', 'horiz', 'center');
end
function [sma, S] = PrepareStateMachine(S, TrialTypes, currentTrial, ~)
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
end
if S.context(currentTrial) == 1 % Reinforced context 
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', 2, ...
        'StateChangeConditions', {'Tup', 'Stimulus', 'Port1In', 'Stop'}, ... % if no action, stimulus activated; if action, stop period
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'Stop', ... % Stop period to assure no activity during Pre-Trial
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'Dead'}, ...
        'OutputActions', {'SoftCode', Stimulus, 'ValveState', 2});  
    sma = AddState(sma, 'Name', 'Dead', ... % 100ms dead period
        'Timer', .1, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', 2, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 1}); 
    sma = AddState(sma, 'Name', 'Drinking', ... % 4 seconds for drinking
        'Timer', 4,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {}); 
    sma = AddState(sma, 'Name', 'Miss', ... % 2 second time out for miss
        'Timer', 2,...
        'StateChangeConditions', {'Tup', 'exit'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', 2,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'Punish', ... % 7 second punish state
        'Timer', 7,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out 
else % Probe context
    sma = AddState(sma, 'Name', 'PreTrial', ... % Pre trial period ensuring no activity for 2 seconds
        'Timer', 2, ...
        'StateChangeConditions', {'Tup', 'Stimulus', 'Port1In', 'Stop'}, ... % if no action, stimulus activated; if action, stop period
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Stop', ... % Stop period to assure no activity during Pre-Trial
        'Timer', 0, ...
        'StateChangeConditions', {'Tup', 'PreTrial'}, ... % returns to PreTrial
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Stimulus', ... % Tone
        'Timer', S.GUI.SoundDuration, ...
        'StateChangeConditions', {'Tup', 'Dead'}, ...
        'OutputActions', {'SoftCode', Stimulus, 'ValveState', 2}); % Keeps lick tube out 
    sma = AddState(sma, 'Name', 'Dead', ... % 100ms dead period
        'Timer', .1, ...
        'StateChangeConditions', {'Tup', 'WaitForLick'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'WaitForLick', ... % Response period
        'Timer', 2, ... 
        'StateChangeConditions', {'Port1In', INResponse, 'Tup', NOResponse},... % If lick, then open valve to reward. If no lick, miss period
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'OpenValve', ... % Open valve for reward
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Drinking', ... % 4 seconds for drinking
        'Timer', 4,...
        'StateChangeConditions', {'Tup', 'ITI'},...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Miss', ... % 2 second time out for miss
        'Timer', 2,...
        'StateChangeConditions', {'Tup', 'exit'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'CorrectReject', ... % 2 second correct reject state
        'Timer', 2,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'Punish', ... % 7 second punish state
        'Timer', 7,...
        'StateChangeConditions', {'Tup', 'ITI'}, ...
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {'ValveState', 2}); % Keeps lick tube out
end

function UpdateOutcomePlot(TrialTypes, Data)
global BpodSystem
Outcomes = zeros(1,Data.nTrials);% Creates a vector for each completed trial, listing outcomes
for x = 1:Data.nTrials
    if ~isnan(Data.RawEvents.Trial{x}.States.OpenValve)
        Outcomes(x) = -1; % green circle for hits
    elseif TrialTypes(x) == 1 % go
        Outcomes(x) = 1; % green x for Miss
    elseif TrialTypes(x) == 2 & ~isnan(Data.RawEvents.Trial{x}.States.CorrectReject) % If no-go and correct reject
        Outcomes(x) = 2; % unfilled green circle for Correct Rejects
    elseif TrialTypes(x) == 2 % no go.
        Outcomes(x) = 0; % red X for false alarm
    end
end
TrialTypeOutcomePlot(BpodSystem.GUIHandles.OutcomePlot,'update',Data.nTrials+1,TrialTypes,Outcomes);