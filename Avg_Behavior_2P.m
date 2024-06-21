
mouses = {'cd016','cd017','cd036','cd037','cd042','cd044'};
disks = [4;4;5;5;4;4];

nMice = length(mouses);

reinforced = cell(nMice,2);
probe = cell(nMice,2);
for i=1:nMice
    results = Behavior_2P_Auditory_GNG(mouses{i},disks(i));
    close all;
    reinforced(i,:) = results(1,:);
    probe(i,:) = results(2,:);
end

%% Plot averages curves

nDmax = max(cell2mat(reinforced(:,1)));
maxYlim = 5;

figure;
subplot(2,3,1); hold on;
for i=1:nMice
    plot(reinforced{i,1},reinforced{i,2},'k');
end

allreinforced = nan(nDmax,nMice);
for i=1:nMice
    allreinforced(reinforced{i,1},i) = reinforced{i,2};
end
avgreinforced = nanmean(allreinforced,2);
subplot(2,3,1); hold on;
plot(1:nDmax,avgreinforced,'k-','linewidth',2);
ylim([-2 maxYlim]);

subplot(2,3,2); hold on;
for i=1:nMice
    plot(probe{i,1},probe{i,2},'color',[0.8 0.8 0.8]);
end

allprobe = nan(nDmax,nMice);
for i=1:nMice
    allprobe(probe{i,1},i) = probe{i,2};
end
avgprobe = nanmean(allprobe,2);
subplot(2,3,2); hold on;
plot(1:nDmax,avgprobe,'color',[0.8 0.8 0.8],'linewidth',2);
ylim([-2 maxYlim]);

subplot(2,3,3); hold on;
plot(1:nDmax,avgprobe,'color',[0.8 0.8 0.8],'linewidth',2);
plot(1:nDmax,avgreinforced,'k','linewidth',2);
ylim([-2 maxYlim]);

allprobecorrrected = nan(nDmax,nMice);
subplot(2,3,4); hold on;
for i=1:nMice
    [~,wmax] = max(probe{i,2});
    probecorrected = probe{i,2};
    for j=2:length(probecorrected)
        if probecorrected(j)<=max(probecorrected(1:(j-1))) 
            probecorrected(j) = nan; % remove worst than any previous day
        end
    end    
    probecorrected(wmax:end) = probecorrected(wmax); % fix the performance
    allprobecorrrected(probe{i,1},i) = probecorrected;
    ok = ~isnan(probecorrected);
    plot(probe{i,1}(ok),probecorrected(ok),'.-','color',[0.8 0.8 0.8]);
end

subplot(2,3,4); hold on;
plot(1:nDmax,avgprobecorrected,'color',[0.8 0.8 0.8],'linewidth',2);
ylim([-2 maxYlim]);

avgprobecorrected = nanmean(allprobecorrrected,2);
semprobecorrected = nansem(allprobecorrrected')';
semreinforced = nansem(allreinforced')';
subplot(2,3,5); hold on;
plot(1:nDmax,avgprobecorrected,'color',[0.8 0.8 0.8],'linewidth',2);
plot(1:nDmax,avgprobecorrected+semprobecorrected,'--','color',[0.8 0.8 0.8]);
plot(1:nDmax,avgprobecorrected-semprobecorrected,'--','color',[0.8 0.8 0.8]);

plot(1:nDmax,avgreinforced,'k','linewidth',2);
plot(1:nDmax,avgreinforced+semreinforced,'k--');
plot(1:nDmax,avgreinforced-semreinforced,'k--');

ylim([-2 maxYlim-1]);

figure;
for i=1:nMice
    subplot(2,3,i); hold on;
    x = (1:nDmax)';
    ok = ~isnan(allprobe(:,i));
    plot(x(ok),allprobe(ok,i),'.-','color',[0.8 0.8 0.8]);
    plot(x,allreinforced(:,i),'k.-');
    ylim([-1.5 maxYlim]);
    PlotHVLines(1.5,'h','k:');
    title(mouses{i});
end

figure;
for i=1:nMice
    subplot(2,3,i); hold on;
    x = (1:nDmax)';
    ok = ~isnan(allprobecorrrected(:,i));
    plot(x(ok),allprobecorrrected(ok,i),'.-','color',[0.8 0.8 0.8]);
    plot(x,allreinforced(:,i),'k.-');
    ylim([-1.5 maxYlim]);
    PlotHVLines(1.5,'h','k:');
    title(mouses{i});
end

% figure;
% for i=1:nMice
%     subplot(2,3,i);hold on;
%     plot(allprobecorrrected(:,i),'.','markersize',10,'color',[0.8 0.8 0.8]);
% end



% [d pd] = allfitdist(avgprobe);
% if size(d(1,1).Params,2)==2,
%     y = cdf(d(1,1).DistName,avgprobe,d(1,1).Params(1),d(1,1).Params(2));
% elseif size(d(1,1).Params,2)==3,
%     y = cdf(d(1,1).DistName,avgprobe,d(1,1).Params(1),d(1,1).Params(2),d(1,1).Params(3));
% end




% figure;plot(1:nDmax,y,'r','linewidth',2);
% hold on;plot(1:nDmax,avgprobe,'k','linewidth',2);
% 
% [sort_x, idx] = sort(probecorrected);
% subplot(2,1,i);hold on;ylabel('Cumulative fraction');
% xlabel('Error (m)');title(behavNames{i});
% plot(sort_x,y(idx),colors{cond});
%     
%     
% D(1).ParamNames{:}
% pd = fitdist(probecorrected,'burr');
% Y = evalmf(mf,X);
% 
% [d pd] = allfitdist(probecorrected);
% d(1).DistName
% pd{1}
% pd = fitdist(probecorrected,d(1).DistName)
% 
% pd = fitdist(probecorrected,'ExtremeValue');
% xgrid = 1:nDmax;
% pdfEst = pdf(pd,xgrid);
% 
% line(xgrid,pdfEst)
% 
% pd = makedist('ExtremeValue','mu',3.04085,'sigma',0.696998)
% pdf(pd)
