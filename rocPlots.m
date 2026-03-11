function []=rocPlots(mgbTempTestsOnly,reinfcolor,optocolor)

dpF = mgbTempTestsOnly{2,39};
dpT = mgbTempTestsOnly{2,40};
dpC = mgbTempTestsOnly{2,41};

crF = mgbTempTestsOnly{4,39};
crT = mgbTempTestsOnly{4,40};
crC = mgbTempTestsOnly{4,41};

hitFaF = mgbTempTestsOnly{6,39};
hitFaT = mgbTempTestsOnly{6,40};
hitFaC = mgbTempTestsOnly{6,41};

dprimeFig=figure; %hit and false alarm rate scatter by animal
subplot(1,3,1);scatter(hitFaF(2,:),hitFaF(1,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
hold on;
scatter(hitFaF(4,:),hitFaF(3,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
ylim([0 1]); xlim([0 1]); title('Full'); ylabel('Hit Rate'); xlabel('FA Rate'); 

subplot(1,3,2);
scatter(hitFaT(2,:),hitFaT(1,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
hold on;
scatter(hitFaT(4,:),hitFaT(3,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
ylim([0 1]); xlim([0 1]);title('Tone'); ylabel('Hit Rate'); xlabel('FA Rate'); 

subplot(1,3,3);
scatter(hitFaC(2,:),hitFaC(1,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
hold on;
scatter(hitFaC(4,:),hitFaC(3,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
ylim([0 1]); xlim([0 1]);title('Choice'); ylabel('Hit Rate'); xlabel('FA Rate'); 
dprimeFig.Position(3:4)=[525 175];
saveas(gcf,['By Animal_hitFaScatter_Opto']);
saveas(gcf,['By Animal_hitFaScatter_Opto.png']);
saveas(gcf,['By Animal_hitFaScatter_Opto.pdf']);

clear dprimeFig
dprimeFig=figure; %criterion by d' 
subplot(1,3,1);
scatter(crF(1,:),dpF(1,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
hold on;
scatter(crF(2,:),dpF(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title('Full'); xlabel('Criterion'); ylabel('d prime'); xlim([-1 2]);ylim([0 4]);

subplot(1,3,2);
scatter(crT(1,:),dpT(1,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
hold on;
scatter(crT(2,:),dpT(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title('Tone'); xlabel('Criterion'); ylabel('d prime'); xlim([-1 2]);ylim([0 4]);

subplot(1,3,3);
scatter(crC(1,:),dpC(1,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
hold on;
scatter(crC(2,:),dpC(2,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title('Choice'); xlabel('Criterion'); ylabel('d prime'); xlim([-1 2]);ylim([0 4]);

dprimeFig.Position(3:4)=[525 175];
saveas(gcf,['By Animal_dpCritScatter_Opto']);
saveas(gcf,['By Animal_dpCritScatter_Opto.png']);
saveas(gcf,['By Animal_dpCritScatter_Opto.pdf']);


% action rate conservation
rfAction=(hitFaF(2,:)+hitFaF(1,:))/2;
fAction=(hitFaF(3,:)+hitFaF(4,:))/2;
rtAction=(hitFaT(2,:)+hitFaT(1,:))/2;
tAction=(hitFaF(3,:)+hitFaF(4,:))/2;
rcAction=(hitFaC(2,:)+hitFaC(1,:))/2;
cAction=(hitFaC(3,:)+hitFaC(4,:))/2;

rateFig=figure;
subplot(1,2,1);
gg=bar([mean(rfAction) mean(fAction) mean(rtAction) mean(tAction) ...
    mean(rcAction) mean(cAction)]);
gg(1).FaceColor='flat';gg(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor;reinfcolor;optocolor];hold on;
ylim([0 0.75]); ylabel('Action rate'); 
scatter(repmat(gg(1).XEndPoints(1),size(rfAction,2),1), ...
    rfAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(gg(1).XEndPoints(2),size(fAction,2),1), ...
    fAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(gg(1).XEndPoints(3),size(rtAction,2),1), ...
    rtAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(gg(1).XEndPoints(4),size(tAction,2),1), ...
    tAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(gg(1).XEndPoints(5),size(rcAction,2),1), ...
    rcAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(gg(1).XEndPoints(6),size(cAction,2),1), ...
    cAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(mean(rfAction),mean(fAction));
sigstar({[1,2]}, p)
[h,p,ci,stats] = ttest2(mean(rtAction),mean(tAction));
sigstar({[3,4]}, p)
[h,p,ci,stats] = ttest2(mean(rcAction),mean(cAction));
sigstar({[5,6]}, p);
box off;

subplot(1,2,2);
ff=bar([mean(rfAction)-mean(fAction) mean(rtAction)-mean(tAction) ...
    mean(rcAction)-mean(cAction)]); box off;
ff(1).FaceColor='flat';ff(1).CData=[optocolor;optocolor;optocolor];hold on;
ylim([0 0.3]); ylabel('Action Rate, light off - light on'); 
scatter(repmat(gg(1).XEndPoints(1),size(rfAction,2),1), ...
    rfAction-fAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(gg(1).XEndPoints(2),size(rtAction,2),1), ...
    rtAction-tAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(gg(1).XEndPoints(3),size(rtAction,2),1), ...
    rcAction-cAction,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
[h,p,ci,stats] = ttest2(mean(rfAction)-mean(fAction),...
    mean(rtAction)-mean(tAction));
sigstar({[1,2]}, p)
[h,p,ci,stats] = ttest2(mean(rtAction)-mean(tAction),...
     mean(rcAction)-mean(cAction));
sigstar({[2,3]}, p)
[h,p,ci,stats] = ttest2(mean(rfAction)-mean(fAction),...
    mean(rcAction)-mean(cAction));
sigstar({[1,3]}, p)
xticklabels({'Full','Stimulus','Choice'});
rateFig.Position(3:4)=[425 175];
saveas(gcf,['By Animal_rateBar_Opto']);
saveas(gcf,['By Animal_rateBar_Opto.png']);
saveas(gcf,['By Animal_rateBar_Opto.pdf']);

%probe behavior 
pFHit=[];pFFa=[];pTHit=[];pTFa=[];pRHit=[];
pRFa=[];
% find days to compare probe
days{1,6}='ExpertDays';
for oo=2:size(days,1)
    if oo==2
        days{oo,6}=[1,2];
    else
        days{oo,6}=[days{oo,3}(1,1)-1, days{oo,3}(1,1)-2];
    end
end

for pp=2:size(days,1)
    subrates=rates{pp,2};
    probeFHit=subrates(days{pp,4},3);
    probeFFa=subrates(days{pp,4},4);
    probeTHit=subrates(days{pp,5},3);
    probeTFa=subrates(days{pp,5},4);
    probeRHit=subrates(days{pp,6},3);
    probeRFa=subrates(days{pp,6},4);
    pFHit=[pFHit;probeFHit];
    pFFa=[pFFa;probeFFa];
    pTHit=[pTHit;probeTHit];
    pTFa=[pTFa;probeTFa];
    pRHit=[pRHit;probeRHit];
    pRFa=[pRFa;probeRFa];
end

probeFig=figure;
subplot(1,2,1);
qq=bar([mean(pRHit) mean(pRFa) mean(pFHit) mean(pFFa)]);
qq(1).FaceColor='flat';qq(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];hold on;
ylabel('Probe action rate'); 
scatter(repmat(qq(1).XEndPoints(1),size(pRHit,2),1), ...
    pRHit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qq(1).XEndPoints(2),size(pRFa,2),1), ...
    pRFa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qq(1).XEndPoints(3),size(pFHit,2),1), ...
    pFHit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qq(1).XEndPoints(4),size(pFFa,2),1), ...
    pFFa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title('Full Trial'); box off;
xticklabels({'Hit','FA','Hit','FA'});

clear qq
subplot(1,2,2);
qq=bar([mean(pRHit) mean(pRFa) mean(pTHit) mean(pTFa)]);
qq(1).FaceColor='flat';qq(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];hold on;
ylabel('Probe action rate'); 
scatter(repmat(qq(1).XEndPoints(1),size(pRHit,2),1), ...
    pRHit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qq(1).XEndPoints(2),size(pRFa,2),1), ...
    pRFa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
scatter(repmat(qq(1).XEndPoints(3),size(pTHit,2),1), ...
    pTHit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
scatter(repmat(qq(1).XEndPoints(4),size(pTFa,2),1), ...
    pTFa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
title('Stimulus');
xticklabels({'Hit','FA','Hit','FA'});
probeFig.Position(3:4)=[425 175]; box off;
saveas(gcf,['By Animal_probeBar_conditionOpto']);
saveas(gcf,['By Animal_probeBar_conditionOpto.png']);
saveas(gcf,['By Animal_probeBar_conditionOpto.pdf']);
end






