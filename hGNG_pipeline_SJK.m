clear
load('datastructAW_88.mat')
for i=1:length(DataStructure)
    DataStructure(i).Performance=sum(DataStructure(i).pTWcorr);
end
T = struct2table(DataStructure);
sortedT = sortrows(T, 'Performance');
DataStructure = table2struct(sortedT);
%% Action Rates
% for i=1:length(DataStructure)
%     DataStructure(i).stimulipTW=DataStructure(i).Stimuli(find(DataStructure(i).Context==1));
%     DataStructure(i).stimulipFF=DataStructure(i).Stimuli(find(DataStructure(i).Context==2));
% end
for i=1:length(DataStructure)
    for ii=1:length(DataStructure(i).pFFkeys)
        if DataStructure(i).pFFkeys(ii)=='["space"]' && DataStructure(i).pFFcorr(ii)==1 %%hit
            outcome(i,ii)=1;
        elseif DataStructure(i).pFFkeys(ii)=='none' && DataStructure(i).pFFcorr(ii)==0 %%miss
            outcome(i,ii)=2;
        elseif DataStructure(i).pFFkeys(ii)=='none' && DataStructure(i).pFFcorr(ii)==1 %%correct reject
            outcome(i,ii)=3;
        elseif DataStructure(i).pFFkeys(ii)=='["space"]' && DataStructure(i).pFFcorr(ii)==0 %%false alarm space
            outcome(i,ii)=4;
        else
            outcome(i,ii)=5; %%error
        end
    end
    for ii=1:10
        HitRate(i,ii)=sum(outcome(i,[1:4]+(ii-1)*4)==1)/2;
        MissRate(i,ii)=sum(outcome(i,[1:4]+(ii-1)*4)==2)/2;
        CrRate(i,ii)=sum(outcome(i,[1:4]+(ii-1)*4)==3)/2;
        FaRate(i,ii)=sum(outcome(i,[1:4]+(ii-1)*4)==4)/2;
    end
end
HitRate(55,3)=1;
HitMean=nanmean(HitRate);
HitSeom=nanstd(HitRate)/sqrt(size(outcome,2));
MissMean=nanmean(MissRate);
MissSeom=nanstd(MissRate)/sqrt(size(outcome,2));
CrMean=nanmean(CrRate);
CrSeom=nanstd(CrRate)/sqrt(size(outcome,2));
FaMean=nanmean(FaRate);
FaSeom=nanstd(FaRate)/sqrt(size(outcome,2));
GoErrorRate=1-(HitRate+MissRate);
GoErrorMean=nanmean(GoErrorRate);
GoErrorSeom=nanstd(GoErrorRate)/sqrt(size(outcome,2));
NoGoErrorRate=1-(CrRate+FaRate);
NoGoErrorMean=nanmean(NoGoErrorRate);
NoGoErrorSeom=nanstd(NoGoErrorRate)/sqrt(size(outcome,2));
colorsmat2=[46, 196, 182; 151, 178, 105; 255, 159, 28; 243, 108, 49; 230, 57, 70; 29, 53, 87; 230, 57, 70; 29, 53, 87]./255;
figure
% shadedErrorBar([1:10],HitMean,HitSeom,'lineProps',{'color', colorsmat2(1,:),'LineWidth',2});
% shadedErrorBar([1:10],FaMean,FaSeom,'lineProps',{'color', colorsmat2(3,:),'LineWidth',2});
% shadedErrorBar([1:10],GoErrorMean,GoErrorSeom,'lineProps',{'color', colorsmat2(2,:),'LineWidth',2});
% shadedErrorBar([1:10],NoGoErrorMean,NoGoErrorSeom,'lineProps',{'color', colorsmat2(4,:),'LineWidth',2});
legend('Hit','FA','TargetError','FoilError')
ylim([0,1])
figure
imagesc(HitRate);
title('HitRate')
figure
imagesc(CrRate);
title('CrRate')
figure
imagesc(GoErrorRate);
title('TargetErrorRate')
figure
imagesc(NoGoErrorRate);
title('FoilErrorRate')
figure
imagesc((HitRate+CrRate)/2);
title('Performance')
%% k means
AvgHitperPerson=mean(HitRate,2);
AvgFAperPerson=mean(FaRate,2);
AvgGoErrorperPerson=mean(GoErrorRate,2);
AvgNoGoErrorperPerson=mean(NoGoErrorRate,2);
% ktest=[AvgHitperPerson,AvgFAperPerson,AvgGoErrorperPerson,AvgNoGoErrorperPerson];
ktest=[HitRate,FaRate,GoErrorRate,NoGoErrorRate];
cluster=4
idx = kmeans(ktest,cluster);
%
% for i=1:cluster
%     figure
%     subplot(2,2,1)
%     imagesc(HitRate(find(idx==i),:)); caxis([0,1])
%     title('HitRate')
%     subplot(2,2,2)
%     imagesc(FaRate(find(idx==i),:)); caxis([0,1])
%     title('FaRate')
%     subplot(2,2,3)
%     imagesc(GoErrorRate(find(idx==i),:)); caxis([0,1])
%     title('TargetErrorRate')
%     subplot(2,2,4)
%     imagesc(NoGoErrorRate(find(idx==i),:)); caxis([0,1])
%     title('FoilErrorRate')
% end

figure(3)
for i=1:cluster
    subplot(2,2,1)
    imagesc(HitRate(find(idx==i),:)); caxis([0,1])
    title('HitRate')
    subplot(2,2,2)
    imagesc(FaRate(find(idx==i),:)); caxis([0,1])
    title('FaRate')
    subplot(2,2,3)
    imagesc(GoErrorRate(find(idx==i),:)); caxis([0,1])
    title('TargetErrorRate')
    subplot(2,2,4)
    imagesc(NoGoErrorRate(find(idx==i),:)); caxis([0,1])
    title('FoilErrorRate')
end
for i=1:length(DataStructure)
    pTWscore(i)=sum(DataStructure(i).pTWcorr);
end
figure
for i=1:cluster
% long = HitRate(:,10)                             % longitude data
% lat = FaRate(:,10)                               % latitude data
% rural = GoErrorRate(:,10)                      % percent rural data
% fatalities = NoGoErrorRate(:,10)                      % fatalities data
subplot(2,3,i)
scatter3(AvgHitperPerson(find(idx==i)),AvgFAperPerson(find(idx==i)),AvgGoErrorperPerson(find(idx==i)),40,AvgNoGoErrorperPerson(find(idx==i)),'filled')    % draw the scatter plot
ax = gca;
ax.XDir = 'reverse';
view(-31,14)
xlabel('Hit Rate')
ylabel('False Alarm Rate')
zlabel('Target Error Rate')
xlim([0,1])
ylim([0,1])
zlim([0,1])
caxis([0,1])
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Foil Error Rate';
end
figure
for i=1:cluster
    paper_dotbarplot_modif(i,((HitRate(find(idx==i),10)+CrRate(find(idx==i),10))./2)',colorsmat2(i,:),5,1)
    hold on
end
%% Damola's Section for X analysis
%% Damola's Section for Y analysis