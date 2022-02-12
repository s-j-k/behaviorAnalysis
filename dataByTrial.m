
%%

% this code takes a data structure and turns it into an array where each
% each row is an animal and each column has different features of the
% animal's performance over day windows
dataCell={};
dataCell(1,1) = {'Mouse'};
dataCell(1,2) = {'Del Hit 1-2 day'}; % What the animal did in response
dataCell(1,3) = {'Hit Rate last day'}; % What the animal did in response
dataCell(1,4) = {'FA Rate first day'}; % Trial number by day
dataCell(1,5) = {'FA rate last day'}; % Trial number by day 
dataCell(1,6) = {'Del dPrime R 1-5 day'};
dataCell(1,7) = {'dPrime last day'};
dataCell(1,8) = {'Day animal 80% R'}; % day as a percentage
dataCell(1,9) = {'Hit rate animal achieves 80% RT'}; % corresponds to column 13
dataCell(1,10) = {'Day animal achieves <30% FA'}; % or trial num
dataCell(1,11) = {'Num Days Animal has >1.5 dPrime'}; % or trial num
dataCell(1,12)= {'Mean lick Onset T day 1'}; 
dataCell(1,13)= {'Mean lick Onset T last day'}; 
% total number of licks data.licksnum
%%
% find the day the animal hits 80% hit rate
ThreshDaysFA=0; 
ThreshDaysH=0; 


counter=1;
for ii=1:length(data) % for each animal
    for jj=1:length(data(ii).expday)-1
        pHit = (data(ii).reinforced{1,jj+1}(1,1)/...
        (data(ii).reinforced{1,jj+1}(1,1)+data(ii).reinforced{1,jj+1}(1,2)));
        pFA = (data(ii).reinforced{1,jj+1}(1,4)/...
        (data(ii).reinforced{1,jj+1}(1,3)+data(ii).reinforced{1,jj+1}(1,4)));
   
        if pFA < 0.3 % find the day the animal hits <30% FA rate
            ThreshDaysFA(counter,1)=ii;
            ThreshDaysFA(counter,2) = pFA;
            ThreshDaysFA(counter,3) = jj;
            ThreshDaysFA(counter,4)=length(data(ii).expday);
            ThreshDaysFA(counter,5)=jj/length(data(ii).expday);
        else
        end
        counter=counter+1;
        
        if pHit > 0.8
            ThreshDaysH(counter,1)=ii;
            ThreshDaysH(counter,2) = pHit; % the experiment days where the animal
                                  % hat a >80 hit rate
            ThreshDaysH(counter,3) = jj; 
            ThreshDaysH(counter,4)=length(data(ii).expday);
            ThreshDaysH(counter,5)=jj/length(data(ii).expday);
        else

        end
        
    end
end


ThreshDaysFA = ThreshDaysFA(any(ThreshDaysFA,2),:);
ThreshDaysH = ThreshDaysH(any(ThreshDaysH,2),:); % take out all rows with only 0 


% find all of the days where the animal has a dprime over 1.5
dPrimeDays=0;
for qq=1:length(data)
    idx = find([data(qq).dprime{1, :}] > 1.5);
    dPrimeDays(qq)=length(idx);
end

%now get just the first day they reach this

%% get the indices of the first row of the ThreshDaysH variable
%this is the first day the animal reaches the 80% Hit rate crtierai

targetval=0;
ii=1;
indices=0; % this is for the hit rate
for ii=1:length(data)
    targetval=targetval + 1;
    indices(end+1)=find(ThreshDaysH(:,1)==targetval,1,'first');
end
indices=indices(2:length(indices));
indices=indices';
hitIndices=indices;


targetval=0;
ii=1;
indices=0; % this is for the FA rate
for ii=1:length(data)
    targetval=targetval + 1;
    try 
        indices(end+1)=find(ThreshDaysFA(:,1)==targetval,1,'first');
    catch
        indices(end+1)=NaN;
    end
end
indices=indices(2:length(indices));
indices=indices';
FAIndices=indices;


%%
%%input all calculated values into the dataCell



for ii =1:length(data)
    for jj=1:length(data(ii).expday)-1
    dataCell(ii+1,1)=cellstr(data(ii).mouse);  
    dataCell(ii+1,2)=num2cell((data(ii).reinforced{1,jj+1}(1,1)/...
        (data(ii).reinforced{1,jj+1}(1,1)+data(ii).reinforced{1,jj+1}(1,2)))-...
        (data(ii).reinforced{1,jj}(1,1)/...
        (data(ii).reinforced{1,jj}(1,1)+data(ii).reinforced{1,jj}(1,2))));
    
    % hit rate last day
    dataCell(ii+1,3)=num2cell((data(ii).reinforced{1,end}(1,1)/...
        (data(ii).reinforced{1,end}(1,1)+data(ii).reinforced{1,end}(1,2))));
    
    % FA rate first day
    dataCell(ii+1,4)=num2cell((data(ii).reinforced{1,1}(1,4)/...
        (data(ii).reinforced{1,1}(1,3)+data(ii).reinforced{1,1}(1,4))));
    
    % FA rate last day
    dataCell(ii+1,5)=num2cell((data(ii).reinforced{1,end}(1,4)/...
        (data(ii).reinforced{1,end}(1,3)+data(ii).reinforced{1,end}(1,4))));
    
    % del dPrime R 1-5 day
    % (data(ii).dprime{1,expday}(1,1))
    dataCell(ii+1,6)=num2cell((data(ii).dprime{1,5}(1,1))-(data(ii).dprime{1,1}(1,1)));
    
    % del dPrime R last day
    dataCell(ii+1,7)=num2cell((data(ii).dprime{1,end}(1,1)));
    dataCell(ii+1,12)=num2cell(nanmean(data(ii).firstlickonset{1,1}.RT));
    dataCell(ii+1,13)=num2cell(nanmean(data(ii).firstlickonset{1,length(data(ii).expday)}.RT));
    
    end
end

ii=1;
for oo=1:length(hitIndices)
    dataCell(oo+1,8)=num2cell(ThreshDaysH(hitIndices(oo),5));
    dataCell(oo+1,9)=num2cell(ThreshDaysH(hitIndices(oo),2));
    if isnan(FAIndices(oo))
        dataCell(oo+1,10)=num2cell(NaN);
    else
        dataCell(oo+1,10)=num2cell(ThreshDaysFA(FAIndices(oo),5));
    end
    
    dataCell(oo+1,11)=num2cell(dPrimeDays(oo));

end


%%
% get water restricted mice

%CITRIC ACID VS WATER RESTRICTION MICE COHORT 2 (after covid):
% CA={'SM039', 'SM040', 'SM041', 'SM042', 'SM043', 'SM044', 'SM045', 'SM046'};
% dob={'04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020','04/28/2020'};
% genotype={'WR','WR','CA','CA','WR','WR','CA', 'CA'};

%CITRIC ACID VS WATER RESTRICTION MICE COHORT 3 (after covid):
% CA2={'SM049', 'SM050', 'SM051', 'SM052', 'SM053', 'SM054', 'SM055', 'SM056'};
% dob={'06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020','06/30/2020'};
% genotype={'CA','CA','CA','CA','WR','WR','WR','WR'};
% dir('T:\su\DATA\behaviorData')
% load('dataCell.mat')
waterDataCell={};
waterDataCell=dataCell(1:3,:);
waterDataCell(end+1:end+2,:)=dataCell(6:7,:);
waterDataCell(end+1:end+4,:)=dataCell(14:17,:);

dataArray=cell2mat(waterDataCell(2:end,2:end));

dataArrayFull=cell2mat(dataCell(2:end,2:end));


%%
% now save this to a .csv
% writecell(dataCell, 'dataCell.mat')
dataArrayWR=waterDataCell(:,:);
% dataArrayV=dataCell(1:end,:);
writecell(dataArrayWR, 'dataMouseWR.csv')
% save('dataMouse.csv','dataArrayV')




%%

    figure;
for ii=1:10
    subplot(5,2,ii); hold on;
    scatter(size(dataArray,1),dataArray(1:end,ii));hold off;
end


%%
figure;
scatter3(dataArray(:,7),dataArray(:,9),dataArray(:,10))
%%
% look at d' over time

figure(1);hold on;
for uu=1:10
    subplot(1,10,uu);hold on;
    for ww=1:length(data(uu).expday)
        plot(ww,(data(uu).dprime{1,ww}(1,1)));
    end
end
%%
[idx,C] = kmeans(dataArray,3);
% Plot the clusters and the cluster centroids.

figure
gscatter(dataArray(:,1),dataArray(:,2),idx,'bgm')
hold on
plot(C(:,1),C(:,2),'kx')
legend('Cluster 1','Cluster 2','Cluster 3','Cluster Centroid')

%%
idx=0;c=0;
    figure;
    hold on
    
for tt=2:10
    [clust,C] = kmeans(dataArrayFull,tt);
    subplot(1,10,tt-1)
    silhouette(dataArrayFull,clust)
end
% for tt=2:7
%     [clust,C] = kmeans(dataArray,tt);
% 
%     subplot(1,10,tt-1)
%     silhouette(dataArray,clust)
% end
    
%     plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
%     hold on
%     plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
%     plot(C(:,1),C(:,2),'kx',...
%          'MarkerSize',15,'LineWidth',3) 
%     legend('Cluster 1','Cluster 2','Centroids',...
%            'Location','NW')
%     title 'Cluster Assignments and Centroids'
%     hold off
    
%%
[idx,C] = kmeans(dataArrayFull, 3); % 3 number of classes
figure;
% Cc = [0 .9 .75; 1 0 0; 0 0.4 0.4];
scatter3(dataArrayFull(:,9), dataArrayFull(:,3), dataArrayFull(:,7), 15, idx, 'filled');
zlabel('Day animal achieves 80 percent R correct');
ylabel('FA Rate first day');
xlabel('Day animal achieves <30% FA');
%%
figure;
[idxSemi,C] = kmeans(dataArray, 3); % 3 number of classes
scatter3(dataArray(:,7), dataArray(:,3), dataArray(:,9), 15,idxSemi, 'filled');
% scatter3(C(:,1),C(:,2),C(:,3),'kx')
xlabel('Day animal achieves 80 percent R correct');
ylabel('FA Rate first day');
zlabel('Day animal achieves <30% FA');
%%

% are the mice using one strategy initially to build their model of the
% task?
% after they have made a preliminary model, can they actively use the model
% they are building to make decisions?
% or, is there a period of consolidation that is required for the model to
% be updated before it is used?


[coeff,score,latent,tsquared,explained]=pca(dataArray);
figure(1);
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')

%%
[coeff,score,latent,tsquared,explained]=pca(dataArrayFull);
figure(2);
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
%%
figure;
    [idx,C]= kmeans(dataArrayFull,3);
    gscatter(dataArrayFull(:,1),dataArrayFull(:,2),idx,'bgm')
    hold on
    plot(C(:,1),C(:,2),'kx')
    legend('Cluster 1','Cluster 2','Cluster 3','Cluster Centroid')
    
    %%
    % get your predictors
    % then your response variable
    X = dataArray(:,1:11);
    Y = dataArray(:,12); % day animal achieves <30% FA 
    ncomp=7;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp);
    %%
    plot(1:ncomp,cumsum(100*PCTVAR(2,:)),'-bo');
    xlabel('Number of PLS components');
    ylabel('Percent Variance Explained in y');
    
    %%
    figure;
    yfit = [ones(size(X,1),1) X]*BETA;
    residuals = Y - yfit;
    stem(residuals)
    xlabel('Observations');
    ylabel('Residuals');
    
    %%
    % calculate the normalized PLS weights
    W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
    
    p = size(XL,1);
sumSq = sum(XS.^2,1).*sum(YL.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));

indVIP = find(vipScore >= 1);
figure;
% plot the VIP scores
% 3 and 7 predictor variabels
scatter(1:length(vipScore),vipScore,'x')
hold on
scatter(indVIP,vipScore(indVIP),'rx')
plot([1 length(vipScore)],[1 1],'--k')
hold off
axis tight
xlabel('Predictor Variables')
ylabel('VIP Scores')



%%


 for di=1:length(data(ii).expday)
        dprimeR{ii}(di,:)=data(ii).dprime{1,di}(1);
        
        dprimeP{ii}(di,:)=data(ii).dprime{1,di}(2);
 end
 %% plot d; over days
figure;
tt=length(dprimeR(WTtouseyoung));
for kk=1:length(dprimeR(WTtouseyoung))
    subplot(tt,1,kk);
    plot(dprimeR{1,WTtouseyoung(kk)});
end
%% plot d' over days on the same plot
figure; hold on;
offset=10;
for kk = 1:tt
    plot(dprimeR{1,WTtouseyoung(kk)}+kk*offset - offset)
end

%% plot FA over days
figure;
for kk=1:length(FalseAlarmsR(WTtouseyoung))
    plot(FalseAlarmsR{1,WTtouseyoung(kk)});hold on;
end


%%
figure;
for kk=1:length(dprimeR(WTtouseyoung))
    plot(dprimeP{1,WTtouseyoung(kk)});hold on;
end

%%
for ii=1:length(data)
    for di=1:length(data(ii).expday)
        ProbeTrials=find(data(ii).trialhist{1,di}>=3);
        ProbeTTargets=find(data(ii).trialhist{1,di}==3);
        ProbeTFoils=find(data(ii).trialhist{1,di}==4);
        Foils=find(data(ii).trialhist{1,di}==2);
        FalseAlarmsR{ii}(di, :)=data(ii).reinforced{1,di}(4)/length(Foils);
        FalseAlarmsP{ii}(di, :)=data(ii).probe{1,di}(4)/length(ProbeTFoils);
    end
end



%%

% this code takes a data structure and turns it into an array where each
% row is a trial

% first, make a cell array where every row is a trial
dataCell={};
dataCell(1,1) = {'Mouse'};
dataCell(1,2) = {'SessionDay'};
dataCell(1,3) = {'Trial Response'}; % What the animal did in response
%        ReinfHit=length(intersect(ReinfAll, find(TrialResp==1)));
%        ReinfMiss=length(intersect(ReinfAll, find(TrialResp==2)));
%        ReinfCR=length(intersect(ReinfAll, find(TrialResp==3)));
%        ReinfFA=length(intersect(ReinfAll, find(TrialResp==4)));

dataCell(1,4) = {'Trial Number'}; % Trial number by day
dataCell(1,5) = {'Trial Type Hist'}; % What type of trial it was 
% 1 is Target
% 2 is Foil
% 3 is Probe Target
% 4 is probe foil

dataCell(1,6) = {'Correct'};
% dataCell(1,7) = {'Probe'};
% dataCell(1,8) = {'Trial History'};
% dataCell(1,9) = {'First Lick Onset Reinf'};
% dataCell(1,10) = {'First Lick Onset'};
% total number of licks data.licksnum

% dataCell(1,11) = {'Num Licks per Trial'}; % in data structure Data.licksall
% dataCell(1,12) = {'Lick frequency'};

%% make a matrix of all the trials
lastNumTrials = 0;
for ii = 1:length(data)  
    if ii==1
        lastNumTrials = 2;
        for numDay=1:length(data(ii).expday)
        
        numTrials = data(ii).ntrials{1,numDay};
        dataCell(lastNumTrials:lastNumTrials+numTrials-1,1)=cellstr(data(ii).mouse);     
        dataCell(lastNumTrials:lastNumTrials+numTrials-1,2)=num2cell(data(ii).expday{1,numDay});             
        dataCell(lastNumTrials:lastNumTrials+numTrials-1,3)=num2cell(data(ii).trialresp{1,numDay})';
        dataCell(lastNumTrials:lastNumTrials+numTrials-1,4)=num2cell(1:data(ii).ntrials{1,numDay})';
        dataCell(lastNumTrials:lastNumTrials+numTrials-1,5)=num2cell(data(ii).trialhist{1,numDay})';
        
        lastNumTrials=lastNumTrials+numTrials;
        end
    else
        for numDay=1:length(data(ii).expday)
            % for variables stored as vectors of each day
            numTrials = data(ii).ntrials{1,numDay};
            dataCell(lastNumTrials:lastNumTrials+numTrials-1,1)=cellstr(data(ii).mouse);     
            dataCell(lastNumTrials:lastNumTrials+numTrials-1,2)=num2cell(data(ii).expday{1,numDay});             
            dataCell(lastNumTrials:lastNumTrials+numTrials-1,3)=num2cell(data(ii).trialresp{1,numDay})';
            dataCell(lastNumTrials:lastNumTrials+numTrials-1,4)=num2cell(1:data(ii).ntrials{1,numDay})';
            dataCell(lastNumTrials:lastNumTrials+numTrials-1,5)=num2cell(data(ii).trialhist{1,numDay})';
            lastNumTrials=lastNumTrials+numTrials;
        end
        
        
%         %for variables stored as vectors of each trial
%             if isnan(num2cell(data(ii).licksall{1,numDay}))
%                 dataCell(lastNumTrials:lastNumTrials+numTrials-1,12)=10;  %licks all also has NaNs. 
%             else
%                 dataCell(lastNumTrials:lastNumTrials+numTrials-1,12)=num2cell(length(data(ii).licksall{1,1}{1,numDay})/max(data(ii).licksall{1,1}{1,numDay}));
%             end 
%             
% %          need to match certain trial types with the corresponding entry of the row        
% %         dataCell(lastNumTrials:lastNumTrials+numTrials-1,9)=num2cell(1:data(ii).firstlickonset{1,numDay}.RT)';
% %           Problem with the above line is that it has NaNs; this is the first lick of each trial type                   
%         
        
    end
end


for jj = 2:length(dataCell)
    
    if mod(cell2mat(dataCell(jj,3)),2)==1
        dataCell(jj,6)=num2cell(1); % animal got it right
    else
        dataCell(jj,6)=num2cell(0); 
    end
end

%% Borrowed from Sharlen's code
maxsizemat=500;
winincrease3=100;
maxdpsli100R=NaN(length(data),maxsizemat);
for ii=1:length(data)
    BigMat{ii}(1,:)=horzcat(data(ii).trialhist{:});%creates matrix with all info per trials not days for type of trial
    BigMat{ii}(2,:)=horzcat(data(ii).trialresp{:});%creates matrix with all info per trials not days for response
    
    
end
ii=1;
for ii=1:length(data)
    clear ReinforcedTrials ProbeTrials TtrialHistR TtrialHistP TrialResponseR TrialResponseP
    ReinforcedTrials=find(BigMat{ii}(1,:)<=2);
    ProbeTrials=find(BigMat{ii}(1,:)>2);
    TtrialHistR=BigMat{ii}(1,ReinforcedTrials);
    TrialResponseR=BigMat{ii}(2,ReinforcedTrials);
    TtrialHistP=BigMat{ii}(1,ProbeTrials);
    TrialResponseP=BigMat{ii}(2,ProbeTrials);
    sliwi=0; %for reinforced
      for sli=1:((max(ReinforcedTrials)-winincrease3)/winincrease3)
          til=find(ReinforcedTrials>sliwi & ReinforcedTrials<=sliwi+winincrease3);
          tillen{ii}(sli,:)=length(til);
                    Rhitssli=sum(TrialResponseR(til)==1);
            Rmisssli=sum(TrialResponseR(til)==2);
            Rcrsli=sum(TrialResponseR(til)==3);
            Rfasli=sum(TrialResponseR(til)==4);
          [dp_sli100R(ii,sli), c_sli100R(ii,sli)]=dprime_criterion_calc(Rhitssli, Rmisssli, Rcrsli, Rfasli);

          [maxdpsli100R(ii,sli) maxcsli100R(ii,sli)]=dprime_criterion_calc(sum(TtrialHistR(til)==1), 0, sum(TtrialHistR(til)==2), 0);
      end
end

%%

maxday=300;
figure
hold all
for ii=wtmice
plot([1:maxday]*100,(dp_sli100R(ii,1:maxday)*100)./maxdpsli100R(ii,1:maxday), 'k-', 'linewidth', 2)
ylabel('Normalized performance')
xlabel('Trials')
text([200],[80], 'saline', 'Color', 'k','FontSize', 12, 'fontname','arial')
text([200],[70], 'dsp4', 'Color', plotscolor,'FontSize', 12, 'fontname','arial')
end 
for ii=mutmice
plot([1:maxday]*100,(dp_sli100R(ii,1:maxday)*100)./maxdpsli100R(ii,1:maxday), '-', 'color', plotscolor, 'linewidth', 2)
end 

for ii=1:length(mouse)
    maxdpsli100Rnorm(ii,:)=max((dp_sli100R(ii,1:maxday)*100)./maxdpsli100R(ii,1:maxday));
    maxdpsli100Pnorm(ii,:)=max((dp_slitouseP(ii,1:maxday)*100)./3.2897);
end



%%
% fit

% see
FitMaking1110
%%


