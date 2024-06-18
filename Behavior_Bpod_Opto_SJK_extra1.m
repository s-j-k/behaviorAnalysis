%% stuff from Bpod_Opto analysis code I don't understand

%% trial blocs
% %% Trial blocs
% sizeBloc = 20;
% nbTrialsPerDay = nan(nFiles,nSubj);
% for i=1:nSubj
%     nbTrialsPerDay(:,i) = Accumulate(MAT{i}(:,SESS));
% end
% 
% for i=1:nFiles
%     for j=1:nSubj
%         trialBloc() = repelem([1:floor(nbTrialsPerDay(i,j)/sizeBloc)],sizeBloc*ones(floor(nbTrialsPerDay(i,j)/sizeBloc),1));
%         trialBloc(floor(nbTrialsPerDay(i,j)/sizeBloc)*sizeBloc+1:nbTrialsPerDay(i,j)) = nan;
%     end
% end



%%

%     fig = figure; hold on;
%     subplot(2,3,1); hold on; % HIT rate, test vs control
%     title('OPTO trials, TEST vs CTL');
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [test_subjrates_hit(:);ctl_subjrates_hit(:)];
%     barwitherrn(data,groups);hold on;
%     g1 = linspace(0.75,1.25,nSubj/2*nFiles);
%     g2 = linspace(1.75,2.25,nSubj/2*nFiles);
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('HIT rate');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,2); hold on; % FA rate, test vs control
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [test_subjrates_fa(:);ctl_subjrates_fa(:)];
%     barwitherrn(data,groups);hold on;
%     g1 = linspace(0.75,1.25,nSubj/2*nFiles);
%     g2 = linspace(1.75,2.25,nSubj/2*nFiles);
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('FA rate');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,3); hold on; % correct
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [correct_test(:);correct_ctl(:)];
%     barwitherrn(data,groups);hold on;
%     g1 = linspace(0.75,1.25,nSubj/2*nFiles);
%     g2 = linspace(1.75,2.25,nSubj/2*nFiles);
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('Portion correct');xlabel('Days');
%     ylim([0 1]);
% 
%     subplot(2,3,4); hold on; % HIT lick latency, test vs control
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mo_hit = mo_hit(ctl,:);
%     test_mo_hit = mo_hit(test,:);
%     data = [test_mo_hit(:);ctl_mo_hit(:)];
%     barwitherrn(data,groups);hold on;
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('HIT lick latency (s)');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,5); hold on; % FA lick latency, test vs control
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mo_fa = mo_fa(ctl,:);
%     test_mo_fa = mo_fa(test,:);
%     data = [test_mo_fa(:);ctl_mo_fa(:)];
%     barwitherrn(data,groups);hold on;
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('FA lick latency (s)');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,6); hold on; % d primes
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mr_hit = odprime(ctl,:);
%     test_mr_hit = odprime(test,:);
%     data = [test_mr_hit(:);ctl_mr_hit(:)];
%     barwitherrn(data,groups);hold on;
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel("d'");xlabel('Days');
%     ylim([-1 4]);
    
%     figure; hold on; % d primes
%     groups = [ones(nSubj/2,1);2*ones(nSubj/2,1)];
%     ctl_mr_hit = pdprime(ctl,4);
%     test_mr_hit = pdprime(test,4);
%     data = [test_mr_hit(:);ctl_mr_hit(:)];
%     barwitherrn(data,groups);hold on;
%     g1 = linspace(0.75,1.25,nSubj/2);
%     g2 = linspace(1.75,2.25,nSubj/2);
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel("d'");xlabel('Days');
%     ylim([-1 3]);
    
%%
%%%%%%%%%%


%     fig = figure; hold on;
%     subplot(2,3,1); hold on; % HIT rate, test vs control
%     title('PROBE, TEST vs CTL');
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [test_subjrates_hit(:);ctl_subjrates_hit(:)];
%     barwitherrn(data,groups);hold on;
%     g1 = linspace(0.75,1.25,nSubj/2*nFiles);
%     g2 = linspace(1.75,2.25,nSubj/2*nFiles);
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('HIT rate');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,2); hold on; % FA rate, test vs control
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [test_subjrates_fa(:);ctl_subjrates_fa(:)];
%     barwitherrn(data,groups);hold on;
%     g1 = linspace(0.75,1.25,nSubj/2*nFiles);
%     g2 = linspace(1.75,2.25,nSubj/2*nFiles);
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('FA rate');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,3); hold on; % correct
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [correct_test(:);correct_ctl(:)];
%     barwitherrn(data,groups);hold on;
%     g1 = linspace(0.75,1.25,nSubj/2*nFiles);
%     g2 = linspace(1.75,2.25,nSubj/2*nFiles);
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('Portion correct');xlabel('Days');
%     ylim([0 1]);
% 
%     subplot(2,3,4); hold on; % HIT lick latency, test vs control
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mp_hit = mp_hit(ctl,:);
%     test_mp_hit = mp_hit(test,:);
%     data = [test_mp_hit(:);ctl_mp_hit(:)];
%     barwitherrn(data,groups);hold on;
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('HIT lick latency (s)');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,5); hold on; % FA lick latency, test vs control
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mp_fa = mp_fa(ctl,:);
%     test_mp_fa = mp_fa(test,:);
%     data = [test_mp_fa(:);ctl_mp_fa(:)];
%     barwitherrn(data,groups);hold on;
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel('FA lick latency (s)');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,6); hold on; % d primes
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     ctl_mr_hit = pdprime(ctl,:);
%     test_mr_hit = pdprime(test,:);
%     data = [test_mr_hit(:);ctl_mr_hit(:)];
%     barwitherrn(data,groups);hold on;
%     plot(g1,data(groups==1),'o');
%     plot(g2,data(groups==2),'o');
%     set(gca,'xtick', [1 2],'xticklabel',{'TEST' 'CTL'});
%     ylabel("d'");xlabel('Days');
%     ylim([-1 4]);

%%

%     fig = figure; hold on;
%     subplot(2,3,1); hold on; % HIT rate, test vs control
%     title('REINF vs OPTO trials, TEST');
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [test_subjrates_hit_r(:);test_subjrates_hit_o(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = test_subjrates_hit_r(i,:);
%         g2 = test_subjrates_hit_o(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     ylabel('HIT rate');xlabel('Days');
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylim([0 1.5]);
% 
%     subplot(2,3,2); hold on; % FA rate 
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [test_subjrates_fa_r(:);test_subjrates_fa_o(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = test_subjrates_fa_r(i,:);
%         g2 = test_subjrates_fa_o(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     ylabel('FA rate');xlabel('Days');
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylabel('FA rate');xlabel('Days');
%     ylim([0 1.5]);
% 
%     subplot(2,3,3); hold on; % correct
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     data = [correct_test_r(:);correct_test_o(:)];
%     barwitherrn(data,groups);hold on;
%     for i=1:nSubj/2
%         g1 = correct_test_r(i,:);
%         g2 = correct_test_o(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylabel('Portion correct');xlabel('Days');
%     ylim([0 1]);
% 
%     subplot(2,3,4); hold on; % HIT lick latency 
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     test_mr_hit = mr_hit(test,:);
%     test_mo_hit = mo_hit(test,:);
%     data = [test_mr_hit(:);test_mo_hit(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = test_mr_hit(i,:);
%         g2 = test_mo_hit(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     ylabel('HIT lick latency (s)');xlabel('Days');
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylim([0 2]);
% 
%     subplot(2,3,5); hold on; % FA lick latency 
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     test_mr_fa = mr_fa(test,:);
%     test_mo_fa = mo_fa(test,:);
%     data = [test_mr_fa(:);test_mo_fa(:)];
%     barwitherrn(data,groups);hold on;
%     colors = {'b','r','g','b','r','g'};
%     for i=1:nSubj/2
%         g1 = test_mr_fa(i,:);
%         g2 = test_mo_fa(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylabel('FA lick latency (s)');xlabel('Days');
%     ylim([0 2]);
% 
%     subplot(2,3,6); hold on; % d primes
%     groups = [ones(nSubj/2*nFiles,1);2*ones(nSubj/2*nFiles,1)];
%     test_mr_hit = rdprime(test,:);
%     test_mo_hit = odprime(test,:);
%     data = [test_mr_hit(:);test_mo_hit(:)];
%     barwitherrn(data,groups);hold on;
%     for i=1:nSubj/2
%         g1 = test_mr_hit(i,:);
%         g2 = test_mo_hit(i,:);
%         plot(1:2, [g1(:) g2(:)]','o-','color',colors{i});
%     end
%     set(gca,'xtick', [1 2],'xticklabel',{'Light OFF' 'Light ON'});
%     ylabel("d'");xlabel('Days');
%     ylim([-1 4]);