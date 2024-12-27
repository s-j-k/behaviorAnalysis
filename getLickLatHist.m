function getLickLatHist(matrix,nbsubj,subjlist)
    SESS=1;CTXT=2;OUTCOME=4;LICKL=8;
    expertHitLicks=[];
    expertHitLicks1=(matrix(matrix(:,SESS)==1 & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL));
    expertHitLicks=[expertHitLicks;expertHitLicks1];

    figure;
    subplot(2,1,1);
%     histogram(expertHitLicks1,30,'Normalization','cdf');
          histogram(expertHitLicks,100,'Normalization','probability');
    title([subjlist{nbsubj} 'Hit Lick Latency for Reinforced Trials at expert, day 8']);
    ytix = get(gca, 'YTick');
%     xlim([-0.05 0.5]);
    set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
    ylabel('Probability Percentage');xlabel('Time after tone onset (s)');
    subplot(2,1,2)
    histogram(expertHitLicks,100,'Normalization','cdf');
%       histogram(expertHitLicks,100,'Normalization','probability');
%     xlim([-0.05 0.5]);
    title([subjlist{nbsubj} 'Hit Lick Latency for Reinforced Trials at expert, days 8 - 17']);
    ytix = get(gca, 'YTick');ylabel('Probability Percentage');
    set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);xlabel('Time after tone onset (s)');
end