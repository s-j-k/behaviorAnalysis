function lickMatOut=getLickLatHist(matrix,nbsubj,subjlist,days,expRange,optoflag)
    SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
    START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
    hitLicksOff={};pHitLicks={};hitLicksFull={};hitLicksTone={};hitLicksChoice={};
    faLicksOff={};pFaLicks={};faLicksFull={};faLicksTone={};faLicksChoice={};
    lickMatOut={};
    if optoflag==0
        expertLicks=[];
        expertLicks1=(matrix(matrix(:,SESS)==days(jj) & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL));
        expertLicks1=[expertLicks;expertLicks1];
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
        title([subjlist{nbsubj} 'Hit Lick Latency for Reinforced Trials at expert']);
        ytix = get(gca, 'YTick');ylabel('Probability Percentage');
        set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);xlabel('Time after tone onset (s)');
    else
        for jj=1:length(expRange)
            hitLicksOff{jj}=matrix(matrix(:,SESS)==expRange(jj) & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL);
            pHitLicks{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==expRange(jj),LICKL);
            hitLicksFull{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==expRange(jj),LICKL);
            hitLicksTone{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==expRange(jj),LICKL);
            hitLicksChoice{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==expRange(jj),LICKL);
            faLicksOff{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==expRange(jj),LICKL);
            pFaLicks{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==expRange(jj),LICKL);
            faLicksFull{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==expRange(jj),LICKL);
            faLicksTone{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==expRange(jj),LICKL);
            faLicksChoice{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==expRange(jj),LICKL);
        end
        lickMatOut={hitLicksOff, pHitLicks, hitLicksFull, ...
            hitLicksTone, hitLicksChoice, faLicksOff, ...
            pFaLicks, faLicksFull, faLicksTone, faLicksChoice};
    end
    
end