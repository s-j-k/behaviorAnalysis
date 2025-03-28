function LickTraining2
global BpodSystem
global Amount
%% Define parameters
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    S.GUI.RewardAmount = 3.0; % ul
end
Amount = S.GUI.RewardAmount;

%% Define trials
MaxTrials = 1000; % Max Trials
%% Initialize plots
TotalRewardDisplay2('init'); % Initialize reward display
%% Main trial loop
for currentTrial = 1:MaxTrials
    R = GetValveTimes(S.GUI.RewardAmount, 1); ValveTime = R;  % Update reward amountsb
    sma = NewStateMatrix(); % Assemble state matrix
    sma = AddState(sma, 'Name', 'WaitForLick', ... % wait to initiate trial
        'Timer', 2, ... 
        'StateChangeConditions', {'Port1In', 'OpenValve', 'SoftCode1', 'ManualDelivery'},... % lick will lead to reward OR button press will lead to manual delivery
        'OutputActions', {'ValveState', 0}); 
    sma = AddState(sma, 'Name', 'OpenValve', ... % water delivery
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 3}); 
    sma = AddState(sma, 'Name', 'ManualDelivery', ... % water delivery
        'Timer', ValveTime,...
        'StateChangeConditions', {'Tup', 'Drinking'},...
        'OutputActions', {'ValveState', 3});
    sma = AddState(sma, 'Name', 'Drinking', ... % 2 seconds for drinking
        'Timer', 2,...
        'StateChangeConditions', {'Tup', 'ITI', 'SoftCode1', 'ManualDelivery'},...
        'OutputActions', {'ValveState', 2});
    sma = AddState(sma, 'Name', 'ITI', ... % Tup Action
        'Timer', 2,...b
        'StateChangeConditions', {'Tup', 'exit'},... % exits trial
        'OutputActions', {'ValveState', 0});
    SendStateMatrix(sma); % Sends a state machine description to a Bpod 
    RawEvents = RunStateMatrix;
    if ~isempty(fieldnames(RawEvents)) % If trial data was returned
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
        BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
        SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
    end
    HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
    if BpodSystem.Status.BeingUsed == 0
        return
    end
%     if ~isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.OpenValve) % update reward display
    if any(isnan(BpodSystem.Data.RawEvents.Trial{currentTrial}.States.ManualDelivery),2)
        TotalRewardDisplay2('add', (S.GUI.RewardAmount));
        TotalRewardDisplay2('licks');
    end
end
end