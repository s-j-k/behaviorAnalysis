function lickMatOut=getLickLatHist(matrix,nbsubj,subjlist,days,optoflag)
    SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
    START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
    hitLicksOff={};pHitLicks={};hitLicksFull={};hitLicksTone={};hitLicksChoice={};
    faLicksOff={};pFaLicks={};faLicksFull={};faLicksTone={};faLicksChoice={};
    lickMatOut={};
    for jj=1:length(days)
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
            hitLicksOff{jj}=matrix(matrix(:,SESS)==days(jj) & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL);
            pHitLicks{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==days(jj),LICKR);
            hitLicksFull{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==days(jj),LICKR);
            hitLicksTone{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==days(jj),LICKR);
            hitLicksChoice{jj}=matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==days(jj),LICKR);
            faLicksOff{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==days(jj),LICKR);
            pFaLicks{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==days(jj),LICKR);
            faLicksFull{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==days(jj),LICKR);
            faLicksTone{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==days(jj),LICKR);
            faLicksChoice{jj}=matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==days(jj),LICKR);
        end
        lickMatOut={hitLicksOff, pHitLicks, hitLicksFull, ...
            hitLicksTone, hitLicksChoice, faLicksOff, ...
            pFaLicks, faLicksFull, faLicksTone, faLicksChoice};
    end
    
end