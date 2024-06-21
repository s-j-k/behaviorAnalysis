range=[1:274];

rhit = sum(tempMat(range,OUTCOME)==1 & tempMat(range,CTXT)==2) / sum(tempMat(range,CTXT)==2 & tempMat(range,TONE)==1); % hit in reinforced light off
phit = sum(tempMat(range,OUTCOME)==1 & tempMat(range,CTXT)==0) / sum(tempMat(range,CTXT)==0 & tempMat(range,TONE)==1); % hit in probe
ohit = sum(tempMat(range,OUTCOME)==1 & tempMat(range,CTXT)==1) / sum(tempMat(range,CTXT)==1 & tempMat(range,TONE)==1); % hit in opto (i.e. reinforced light on)
rfhit = sum(tempMat(range,OUTCOME)==1 & (tempMat(range,CTXT)==1 | tempMat(range,CTXT)==2 | tempMat(range,CTXT)==5 | tempMat(range,CTXT)==6) / sum((tempMat(range,CTXT)==1 | tempMat(range,CTXT)==2 | tempMat(range,CTXT)==5 | tempMat(range,CTXT)==6) & tempMat(range,TONE)==1)); % hit in reinforced light off + on
phit2_1 = sum(tempMat(:,OUTCOME)==1 & probeIdx1) / sum(probeIdx1 & tempMat(:,TONE)==1); % hit in probe, dividing probe bloc in 2
phit2_2 = sum(tempMat(:,OUTCOME)==1 & probeIdx2) / sum(probeIdx2 & tempMat(:,TONE)==1); % hit in probe, dividing probe bloc in 2
phit2 = [phit2_1 phit2_2];
othit = sum(tempMat(range,OUTCOME)==1 & tempMat(range,CTXT)==5) / sum(tempMat(range,CTXT)==5 & tempMat(range,TONE)==1); % hit in opto tone (i.e. reinforced light on)
ochit = sum(tempMat(range,OUTCOME)==1 & tempMat(range,CTXT)==6) / sum(tempMat(range,CTXT)==6 & tempMat(range,TONE)==1); % hit in opto choice (i.e. reinforced light on)

rfa = sum(tempMat(range,OUTCOME)==3 & tempMat(range,CTXT)==2) / sum(tempMat(range,CTXT)==2 & tempMat(range,TONE)==2);
pfa = sum(tempMat(range,OUTCOME)==3 & tempMat(range,CTXT)==0) / sum(tempMat(range,CTXT)==0 & tempMat(range,TONE)==2);
ofa = sum(tempMat(range,OUTCOME)==3 & tempMat(range,CTXT)==1) / sum(tempMat(range,CTXT)==1 & tempMat(range,TONE)==2);
rffa = sum(tempMat(range,OUTCOME)==3 & (tempMat(range,CTXT)==1 | tempMat(range,CTXT)==2)) / sum((tempMat(range,CTXT)==1 | tempMat(range,CTXT)==2) & tempMat(range,TONE)==2);
pfa2_1 = sum(tempMat(:,OUTCOME)==3 & probeIdx1) / sum(probeIdx1 & tempMat(:,TONE)==2);
pfa2_2 = sum(tempMat(:,OUTCOME)==3 & probeIdx2) / sum(probeIdx2 & tempMat(:,TONE)==2);
pfa2 = [pfa2_1 pfa2_2];
otfa = sum(tempMat(range,OUTCOME)==3 & tempMat(range,CTXT)==5) / sum(tempMat(range,CTXT)==5 & tempMat(range,TONE)==2); % hit in opto tone (i.e. reinforced light on)
ocfa = sum(tempMat(range,OUTCOME)==3 & tempMat(range,CTXT)==6) / sum(tempMat(range,CTXT)==6 & tempMat(range,TONE)==2); % hit in opto tone (i.e. reinforced light on)


cofffa = sum(tempMat(range,OUTCOME)==3 & tempMat(range,CTXT)==3) / sum(tempMat(range,CTXT)==3 & tempMat(range,TONE)==0); % catch light off FA
conffa = sum(tempMat(range,OUTCOME)==3 & tempMat(range,CTXT)==4) / sum(tempMat(range,CTXT)==4 & tempMat(range,TONE)==0); % catch light on FA
vals(2,:) = [rhit rfa phit pfa ohit ofa ...
            rfhit rffa cofffa conffa phit2 pfa2 ...
            othit otfa ochit ocfa];
        
        %%
optocolor=[102/255 178/255 255/255];
reinfcolor= [0.4,0.4,0.4];   
figure;subplot(1,2,1);
sss=bar([vals(2,1) vals(1,15); vals(2,2) vals(2,16)]); %hit tone opto hit fa tone opto FA
sss(1).FaceColor='flat';
sss(2).FaceColor='flat';
sss(1).CData=[reinfcolor;reinfcolor];hold on;
sss(2).CData=[optocolor;optocolor];
scatter(sss(1).XEndPoints(1),size(vals(2,1)),vals(2,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(sss(1).XEndPoints(2),size(vals(2,1)),vals(2,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
    
xticklabels({'hit','fa'});ylabel('Average Action Rate');
title('IC Tone Opto');
subplot(1,2,2);
sss=bar([vals(2,1) vals(1,17); vals(2,2) vals(1,18)]); %hit tone opto hit fa tone opto FA
sss(1).FaceColor='flat';
sss(2).FaceColor='flat';
sss(1).CData=[reinfcolor;reinfcolor];
sss(2).CData=[optocolor;optocolor];hold on;
xticklabels({'hit','fa'});
title('IC Choice Opto');
ylim([0 1]);
