%% List of mice and experiments
mice = {'cd017','cd036','cd019','cd042','cd041','cd044','cd037','zz033'};
exps = [1;1;2;1;2;1;1;2]; % 1=GNG, 2= PASSIVE

%% Plot passsive behavior

gray = [0.8 0.8 0.8];
blocreinf = 100;
sliding = 20;
% rate = 0.4; blocSize = 20;
rate = 0; blocSize = 0;
restriction = [rate blocSize];
plotfig = true;
savefig = false;
pathsave = 'M:\Celine\exci_variables\';
path = 'M:\Celine\exci_variables\';
BLOC = 2;
DAY = 1;
learning = 5; % 4= d'; 5= percent;
PASS = 1;

PASSs = find(exps==2);
nPASS = length(PASSs);
passive_days_perf = nan(15,nPASS);
for m=1:nPASS
    mm = PASSs(m);
    mouse = mice{mm};
        
    behav_day = Behavior_2P_Auditory_GNG_new(mouse,[],'windows','days','fixprobe',true,'passive',true,...
        'plotfig',false);   
    passive_days_perf(1:15,m) = behav_day{PASS}(1:15,learning);    
end

figure;
subplot(2,2,1);hold on;
ntoplot = nanmean(passive_days_perf,2);
plot(1:15,ntoplot)%,'linewidth',2,'k');
plot(nanmean(passive_days_perf,2)+(nansem(passive_days_perf')'),'k');
plot(nanmean(passive_days_perf,2)-(nansem(passive_days_perf')'),'k');
ylim([0.4 1]);
xlim([0 16]);
ylabel('Accuracy');
xlabel('Days');
title(['Performace Passive (n=' num2str(nPASS) ')' ]);
%%
REINF = 1; PROBE = 2;
gray = [0.8 0.8 0.8];
blocreinf = 100;
sliding = 20;
% rate = 0.4; blocSize = 20;
rate = 0; blocSize = 0;
restriction = [rate blocSize];
plotfig = true;
savefig = false;
pathsave = 'M:\Celine\exci_variables\';
path = 'M:\Celine\exci_variables\';
BLOC = 2;
DAY = 1;
learning = 5; % 4= d'; 5= percent;

GNGs = find(exps==1);
GNGs(end) = []; % remove cd037
nGNG = length(GNGs);
allbehav_day_reinf = nan(16,nGNG);
allbehav_day_probe = nan(16,nGNG);
allbehav_trials_reinf = nan(6000,nGNG);
allbehav_trials_probe = nan(6000,nGNG);
allbehav_day_reinf_FA = nan(16,nGNG);
allbehav_day_probe_FA = nan(16,nGNG);
allbehav_day_probe_H = nan(16,nGNG);
allbehav_trial_reinf_FA = nan(6000,nGNG);
for m=1:nGNG
    mm = GNGs(m);
    mouse = mice{mm};
        
    behav_day = Behavior_2P_Auditory_GNG_new(mouse,[],'windows','days',...
        'activityreinf',restriction,'fixprobe',true,'plotfig',false);
    probe_days_dprime = behav_day{PROBE}(:,1);
    dprime_probe_days = behav_day{PROBE}(:,learning);        
    [~,wm] = max(dprime_probe_days(ismember(probe_days_dprime,1:10)));
    dprime_probe_days(wm:end) = dprime_probe_days(wm);
      
    behav_trial = Behavior_2P_Auditory_GNG_new(mouse,[],'windows','blocs','blocreinf',blocreinf,...
        'sliding',sliding,'activityreinf',restriction,'fitted',false,'fixprobe',true,'plotfig',false);
    probe_trials_dprime = behav_trial{PROBE}(:,1);  
    dprime_probe_trials = behav_trial{PROBE}(:,learning);        
    [~,wm] = max(dprime_probe_trials(1:10));
    dprime_probe_trials(wm:end) = dprime_probe_trials(wm);
    
    if plotfig
        fig=figure;
        subplot(2,2,1);hold on;
        plot(behav_day{REINF}(:,1),behav_day{REINF}(:,2),'k');
        plot(behav_day{REINF}(:,1),behav_day{REINF}(:,3),'k--');
        plot(behav_day{PROBE}(:,1),behav_day{PROBE}(:,2),'color',gray);
        plot(behav_day{PROBE}(:,1),behav_day{PROBE}(:,3),'--','color',gray);
        ylim([0 1]);
        xlabel('Days');ylabel('Action rate');
        title([mouse ', action rate']);

        subplot(2,2,2);hold on;
        plot(behav_day{REINF}(:,1),behav_day{REINF}(:,learning),'k');
        nop = isnan(dprime_probe_days);
        plot(behav_day{PROBE}(~nop,1),dprime_probe_days(~nop),'color',gray);
        xlabel('Days');ylabel("d'");
        ylim([0.4 1]);
        title('Performance');
        
        
        results = load([path mouse '-results_nosignals.mat']);
        results = results.results;
        matrix = results{3};
        ctx = results{7};
        if strcmp(mouse,'cd017')
            nop = matrix(:,BLOC)==4; % remove bad trials
        else
            nop = ~ismember(matrix(:,BLOC),matrix(:,BLOC)); % take everything
        end        
        trials = Accumulate(matrix(~nop & ctx(:,1),DAY));

        subplot(2,2,3);hold on;
        plot(behav_trial{REINF}(:,1),behav_trial{REINF}(:,2),'k');
        plot(behav_trial{REINF}(:,1),behav_trial{REINF}(:,3),'k--');
        plot(behav_trial{PROBE}(:,1),behav_trial{PROBE}(:,2),'color',gray);
        plot(behav_trial{PROBE}(:,1),behav_trial{PROBE}(:,3),'--','color',gray);
        ylim([0 1]);
        xlabel('Trials');ylabel('Action rate');
        title([num2str(blocreinf) '-trial bloc, sliding ' num2str(sliding)]);
        PlotHVLines(cumsum(trials),'v','color',gray);
        
        subplot(2,2,4);hold on;
        plot(behav_trial{REINF}(:,1),behav_trial{REINF}(:,learning),'k');
        nop = isnan(dprime_probe_trials);
        plot(behav_trial{PROBE}(~nop,1),dprime_probe_trials(~nop),'color',gray);
        xlabel('Trials');ylabel("d'");
        ylim([0.4 1]);
        title('Performance');       
        PlotHVLines(cumsum(trials),'v','color',gray);
        drawnow;
    end
    
    if savefig
        saveas(fig,[pathsave mouse '-Behavior_reinfrestriction_' num2str(rate) '_' num2str(blocSize) '.pdf']);
        close(fig);
    end
    
    allbehav_day_reinf(behav_day{REINF}(:,1),m) = behav_day{REINF}(:,learning);
    allbehav_day_probe(behav_day{PROBE}(:,1),m) = dprime_probe_days;
    allbehav_trials_reinf(behav_trial{REINF}(:,1),m) = behav_trial{REINF}(:,learning);
    allbehav_trials_probe(behav_trial{PROBE}(:,1),m) = dprime_probe_trials;
    
    allbehav_trial_reinf_FA(behav_trial{REINF}(:,1),m) = behav_trial{REINF}(:,3);
    
    allbehav_day_reinf_FA(behav_day{REINF}(:,1),m) = behav_day{REINF}(:,3);
    allbehav_day_probe_FA(behav_day{PROBE}(:,1),m) = behav_day{PROBE}(:,3);
    allbehav_day_probe_H(behav_day{PROBE}(:,1),m) = behav_day{PROBE}(:,2);
end


% allbehav_day_reinf(16:end,:) = [];
% allbehav_day_probe(16:end,:) = [];

nop = sum(isnan(allbehav_trials_reinf),2)==nGNG;
xtrials_reinf = (1:size(allbehav_trials_reinf,1))';
allbehav_trials_reinf(nop,:) = [];
allbehav_trial_reinf_FA(nop,:) = [];
xtrials_reinf(nop,:) = [];
maxtrial = min(xtrials_reinf(sum(isnan(allbehav_trials_reinf),2)==3));

xtrials_probe = (1:size(allbehav_trials_probe,1))';

%% Find criteria to group mice according to behavior (i.e. TCA alignment)

% Expression
figure;hold on;
subplot(2,2,1);hold on;
plot(allbehav_day_reinf_FA,'.-');
criteria = [0.2;0.4;0.6;0.8];
PlotHVLines(criteria,'h','k:');
xlim([0 size(allbehav_day_reinf_FA,1)+1]);
ylabel('FA');
xlabel('Days');

day_selected = nan(nGNG,1);
for m=1:nGNG
    there = find(allbehav_day_reinf_FA(:,m)<0.2);
    day_selected(m) = there(1);
    plot(day_selected(m),allbehav_day_reinf_FA(day_selected(m),m), 'ko','markersize',10,...
        'linewidth',2);
end

subplot(2,2,2);hold on;
for m=1:nGNG
    plot(allbehav_day_reinf_FA(day_selected(m)-6:day_selected(m)+2,m),'.-');
end
ndays_selected = length(day_selected(m)-6:day_selected(m)+2);
ylabel('FA');
xlabel('Days from expert');
set(gca,'xtick',1:ndays_selected,'xticklabels',num2str((-6:2)'));

subplot(2,2,3);hold on;
for m=1:nGNG
    plot(allbehav_day_reinf(:,m),'.-');
end
ylim([0 4.5]);xlim([0 size(allbehav_day_reinf_FA,1)+1]);
PlotHVLines([0.5;1.5;3.1],'h','k:');

subplot(2,2,4);hold on;
for m=1:nGNG
    plot(allbehav_day_reinf(day_selected(m)-6:day_selected(m)+2,m),'.-');
end
ylabel('d"');
xlabel('Days from expert');
set(gca,'xtick',1:ndays_selected,'xticklabels',num2str((-6:2)'));
legend(mice(GNGs));

% Acquisition
figure;hold on;
subplot(2,2,1);hold on;
plot(allbehav_day_probe_FA,'.-');

subplot(2,2,3);hold on;
plot(allbehav_day_probe_H,'.-');

subplot(2,2,2);hold on;
for m=1:nGNG
    nop = isnan(allbehav_trials_probe(:,m));
    plot(xtrials_probe(~nop),allbehav_trials_probe(~nop,m),'.-');
end
criteria = 1.2;
PlotHVLines(criteria,'h','k:');
xlabel('Trials');
ylabel("d'");

trial_selected = nan(nGNG,1);
for m=1:nGNG
    there = find(allbehav_trials_probe(:,m)>=criteria);
    trial_selected(m) = there(1);
    plot(trial_selected(m),allbehav_trials_probe(trial_selected(m),m), 'ko','markersize',10,...
        'linewidth',2);
end

subplot(2,2,4);hold on;
for m=1:nGNG
    in = InIntervals(xtrials_probe,[trial_selected(m)-500 trial_selected(m)+3000]);
%     n = sum(in);
    nop = isnan(allbehav_trials_probe(:,m));
%     plot(1:n,allbehav_trials_reinf(in,m),'.-');
    plot(xtrials_probe(in & ~nop),allbehav_trials_probe(in & ~nop,m),'.-');
end
PlotHVLines(500,'v','k:');


%% to find criteria: days
% figure;hold on;
% plot(allbehav_day_reinf,'.-');
% criteria = [0.5;2;2.7;3.5];
% PlotHVLines(criteria,'h','k:');
% nPhases = length(criteria);
% day_selected = nan(nPhases,nGNG);
% for m=1:nGNG
%     % <0.5
%     day_selected(1,m) = find(allbehav_day_reinf(:,m)<criteria(1));
%     % 2<x<2.5
%     there = find(allbehav_day_reinf(:,m)>criteria(2) & allbehav_day_reinf(:,m)<criteria(3));
%     day_selected(2,m) = there(1);
%     % 2.5<x<3
%     there = find(allbehav_day_reinf(:,m)>criteria(3) & allbehav_day_reinf(:,m)<criteria(4));
%     day_selected(3,m) = there(1);
%     % x>3.5
%     there = find(allbehav_day_reinf(:,m)>criteria(4));
%     day_selected(4,m) = there(1);
%     
% %     % 3<x<3.5
% %     there = find(allbehav_day_reinf(:,m)>criteria(4) & allbehav_day_reinf(:,m)<criteria(5));
% %     day_selected(4,m) = there(1);
% %     % x>3.5
% %     there = find(allbehav_day_reinf(:,m)>criteria(5));
% %     day_selected(5,m) = there(1);
%     
%     plot(day_selected(:,m),allbehav_day_reinf(day_selected(:,m),m), 'ko','markersize',10,...
%         'linewidth',2);
% end
%% 
nprepost = 20;
figure;hold on;
for m=1:length(mice)
    mouse = mice{m};
    results = load([path mouse '-results_nosignals.mat']);
    results = results.results;
    matrix = results{3};
    ctx = results{7};
    if strcmp(mouse,'cd017')
        nop = matrix(:,BLOC)==4; % remove bad trials
    else
        nop = ~ismember(matrix(:,BLOC),matrix(:,BLOC)); % take everything
    end    
    ok = ~nop & ctx(:,1);
    rmatrix = matrix(ok,:);
    trials = Accumulate(rmatrix(:,DAY));
    trials(trials==0)=[];
    cstrials = cumsum(trials);
    idx = nan(3,nprepost,length(mice));
    for p=1:nPhases
        if p==1
            idx_middle = round(cstrials(day_selected(p,m))/2);
        else
            idx_middle = cstrials(day_selected(p,m)-1)+round(trials(day_selected(p,m))/2);
        end 
        idx_H = rmatrix(ks(ok,H)&rmatrix(:,TRIAL)<idx_middle,TRIAL); 
        idx_H_pre = idx_H(end-nprepost+1:end);
        idx_H = rmatrix(ks(ok,H)&rmatrix(:,TRIAL)>=idx_middle,TRIAL); 
        idx_H_post = idx_H(1:nprepost);                
        idx(H,:,m) = [idx_H_pre;idx_H_post];
        
        idx_FA = rmatrix(ks(ok,FA)&rmatrix(:,TRIAL)<idx_middle,TRIAL); 
        idx_FA_pre = idx_FA(end-nprepost+1:end);
        idx_FA = rmatrix(ks(ok,FA)&rmatrix(:,TRIAL)>=idx_middle,TRIAL); 
        idx_FA_post = idx_FA(1:nprepost);                
        idx(FA,:,m) = [idx_FA_pre;idx_FA_post];
        
        idx_H = rmatrix(ks(ok,H)&rmatrix(:,TRIAL)<idx_middle,TRIAL); 
        idx_H_pre = idx_H(end-nprepost+1:end);
        idx_H = rmatrix(ks(ok,H)&rmatrix(:,TRIAL)>=idx_middle,TRIAL); 
        idx_H_post = idx_H(1:nprepost);                
        idx(H,:,m) = [idx_H_pre;idx_H_post];
    end
        
%         ok2 = ok &

        plot(xtrials_reinf,allbehav_trials_reinf(:,m))

end
%%
fig=figure;
subplot(2,2,1);hold on;
plot(allbehav_day_reinf(:,2),'k');
x = 1:16;
nop = isnan(allbehav_day_probe(:,2));
plot(x(~nop),allbehav_day_probe(~nop,2),'color',gray);
xlabel('Days');ylabel("d'");
ylim([0.4 1]);
title('Performance');

subplot(2,2,2);hold on;
plot(xtrials_reinf,allbehav_trials_reinf(:,2),'k');
nop = isnan(allbehav_trials_probe(:,2));
plot(xtrials_probe(~nop),allbehav_trials_probe(~nop,2),'color',gray);
xlabel('Trials');ylabel("d'");
ylim([0.4 1]);
title('Performance');

%%
savefig = false;
fig =figure;
x = 1:16;
subplot(2,2,1);hold on;
plot(x,nanmean(allbehav_day_reinf,2),'k','linewidth',2);
plot(x,nanmean(allbehav_day_reinf,2)+nansem(allbehav_day_reinf')','k');
plot(x,nanmean(allbehav_day_reinf,2)-nansem(allbehav_day_reinf')','k');
plot(x,nanmean(allbehav_day_probe,2),'color',gray,'linewidth',2);
plot(x,nanmean(allbehav_day_probe,2)+nansem(allbehav_day_probe')','color',gray);
plot(x,nanmean(allbehav_day_probe,2)-nansem(allbehav_day_probe')','color',gray);
xlabel('Days');
% ylim([0 4.5]);title(['Performance (n=' num2str(nGNG) ')']);ylabel("d'");
ylim([0.4 1]);title(['Performance (n=' num2str(nGNG) ')']);ylabel("%");

subplot(2,2,2);hold on;
plot(xtrials_reinf,nanmean(allbehav_trials_reinf,2),'k','linewidth',2);
plot(xtrials_reinf,nanmean(allbehav_trials_reinf,2)+nansem(allbehav_trials_reinf')','k');
plot(xtrials_reinf,nanmean(allbehav_trials_reinf,2)-nansem(allbehav_trials_reinf')','k');
for m=1:nGNG
    nop = isnan(allbehav_trials_probe(:,m));
    plot(xtrials_probe(~nop),allbehav_trials_probe(~nop,m),'.','color',gray,'markersize',10);
%     pause();
end
nop = sum(isnan(allbehav_trials_probe),2)==nGNG;
plot(xtrials_probe(~nop,:), smoothdata(nanmean(allbehav_trials_probe(~nop,:),2),1,'movmean',5),...
     'color',gray,'linewidth',2);
xlim([0 maxtrial]);
xlabel('Trials');
% ylim([0 4.5]);title(['Performance (n=' num2str(nGNG) ')']);ylabel("d'");
ylim([0.4 1]);title(['Performance (n=' num2str(nGNG) ')']);ylabel("%");

subplot(2,2,3);hold on;
for m=1:nGNG
    nop = isnan(allbehav_trials_probe(:,m));
    plot(xtrials_probe(~nop),allbehav_trials_probe(~nop,m),'.','color',gray,'markersize',10);
%     pause();
end

nop = sum(isnan(allbehav_trials_probe),2)==nGNG;
f = @(param,xval) param(1) + ( param(2)-param(1) )./ (   1 + 10.^( ( param(3) - xval ) * param(4) )   );
x = xtrials_probe(~nop,:);
x_vector = min(x):(max(x)-min(x))/100:max(x);        
y = nanmean(allbehav_trials_probe(~nop,:),2);
param = sigm_fit(x,y,[],[],0); 
plot(x,f(param,x),'-','color',gray,'linewidth',2);


avallbehav_trials_probe = nanmean(allbehav_trials_probe,2);
nop = sum(isnan(allbehav_trials_probe),2)==nGNG;
smoothavallbehav_trials_probe = smoothdata(avallbehav_trials_probe,'movmean',500);
figure;hold on;
plot(xtrials_probe,smoothavallbehav_trials_probe,'.');
% fo = fit(x,y,'exp1','StartPoint',[0 0]);
nop = isnan(smoothavallbehav_trials_probe);
fo = fit(xtrials_probe(~nop),smoothavallbehav_trials_probe(~nop),'smoothingspline');
tofit = xtrials_probe(~nop);
plot(tofit,fo(tofit),'linewidth',2);

x = probe_days_dprime;
x_vector=min(x):(max(x)-min(x))/100:max(x);        
y = dprime_probe_days;
param = sigm_fit(x,y,[],[],0); 
plot(x_vector,f(param,x_vector),'-','color',color_probe,'linewidth',2); 
            
            
if savefig
    saveas(fig,[pathsave 'AverageMice-Behavior_reinfrestriction_' num2str(rate) '_' num2str(blocSize) '.pdf']);
    close(fig);
end
