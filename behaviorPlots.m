figure; hold on;
correctidx=[];
for ii = 1:nSubj
    correctidx{ii}=find(rpc(ii,:)>85); %greater than 80% correct in reinf
end

for jj=1:nSubj
    firstHighDay(1,jj)=correctidx{1,jj}(1,1); 
    firstHighDay(2,jj)=rpc(jj,firstHighDay(1,jj));
end

%% plot just percent correct
figure(10);hold on;
avgrpctrace=mean(rpc(:,:));
plot(avgrpctrace,'LineWidth',4)
for pp=1:nSubj
    plot(rpc(pp,:),'Color',[0.7 0.7 0.7]);
end
xlim([1 15]);ylim([50 100]);

%%
figure;
shadedErrorBar(DAYS,nanmean(rpc(:,:)),nansem(rpc(:,:)),{'color','b'},0.5);
xlabel('days');ylabel('percent correct');
xlim([1 15]);ylim([50 100]);
%% 
% plot average d' by days for test and ctl 
shadedErrorBar(DAYS,nanmean(rfdprime(test,:)),nansem(rfdprime(test,:)),{'color','b'},0.5);
shadedErrorBar(DAYS,nanmean(rfdprime(ctl,:)),nansem(rfdprime(ctl,:)),{'color','k'},0.5);

% CDF plot
figure;
cdfplot(firstHighDay(1,:));
ylabel('F(x) proportion of animals');
xlabel('x (first day of >85% correct)')

% normal probability plot
normplot(firstHighDay(1,:));
xlabel('First day of >85% performance');
title('normal probability plot of first day of >85% correct');

% scatter plot of percent correct by days
figure;
scatter(firstHighDay(1,:),firstHighDay(2,:));

