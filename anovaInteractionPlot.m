function anovaInteractionPlot(anovaMat,allDataTestsOnly)
%% Repeated-measures 2 groups (between) x 3 conditions (within)
close all; clc
rng(4);

groups     = ["Test","Ctl"];
conditions = ["Full","Stimulus","Choice"];

nSubPerGroup = 6;     % subjects per group

%% =========================
% Figure 1
%subplot 1 Interaction plot (means ± within-subject CI) for both groups,
% dprime
% =========================

testFullTrial=[allDataTestsOnly{2,39}';allDataTestsOnly{2,40}';allDataTestsOnly{2,41}'];
sigmaOff    = [std(testFullTrial(1:6,1)), std(testFullTrial(7:12,1)), std(testFullTrial(13:18,1))];   % between-subject variability
sigmaOn    = [std(testFullTrial(1:6,2)), std(testFullTrial(7:12,2)), std(testFullTrial(13:18,2))];   % between-subject variability
% True means to create a main effect + interaction 
muA = [mean(testFullTrial(1:6,1)), mean(testFullTrial(7:12,1)), mean(testFullTrial(13:18,1))];    % Group A condition means
muB = [mean(testFullTrial(1:6,2)), mean(testFullTrial(7:12,2)), mean(testFullTrial(13:18,2))];    % Group B condition means (flatter -> interaction)
c2An = multcompare(anovaMat{6,3},"Estimate","row",'Display', 'off'); %full, stimulus, choicce
tbl2AnD = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

reinfcolor= [0.4,0.4,0.4];
optocolor=[102/255 178/255 255/255];
interactionFig=figure;subplot(1,3,1);
errorbar(muA,sigmaOff,'Color',reinfcolor+0.2,'LineWidth',2); hold on;
scatter([1,2,3],muA,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor); hold on;
errorbar(muB,sigmaOn,'Color',optocolor,'LineWidth',2);
scatter([1,2,3],muB,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
xticks([1,2,3]); xticklabels(conditions)
xlim([0.5 3.5]);
ylabel('d prime');box off;
% legend('','Off','','On', 'Location','best')
sigstar({[1,2]}, table2array(tbl2AnD(1,6)))
sigstar({[2,3]}, table2array(tbl2AnD(3,6)))
sigstar({[1,3]}, table2array(tbl2AnD(2,6)))
ylim([-0.5 4]);

%% subplot 2 Interaction plot for both groups hit rate

hitRates=[allDataTestsOnly{6,39}(1,1:6)';allDataTestsOnly{6,40}(1,1:6)';allDataTestsOnly{6,41}(1,1:6)'];
optoHitRates=[allDataTestsOnly{6,39}(3,1:6)';allDataTestsOnly{6,40}(3,1:6)';allDataTestsOnly{6,41}(3,1:6)'];
testFullTrial=horzcat(hitRates, optoHitRates);
sigmaOff    = [std(testFullTrial(1:6,1)), std(testFullTrial(7:12,1)), std(testFullTrial(13:18,1))];   % between-subject variability
sigmaOn    = [std(testFullTrial(1:6,2)), std(testFullTrial(7:12,2)), std(testFullTrial(13:18,2))];   % between-subject variability
% True means to create a main effect + interaction 
muA = [mean(testFullTrial(1:6,1)), mean(testFullTrial(7:12,1)), mean(testFullTrial(13:18,1))];    % Group A condition means
muB = [mean(testFullTrial(1:6,2)), mean(testFullTrial(7:12,2)), mean(testFullTrial(13:18,2))];    % Group B condition means (flatter -> interaction)
c2An = multcompare(anovaMat{4,3},"Estimate","row",'Display', 'off'); %full, stimulus, choicce
tbl2AnD = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
subplot(1,3,2);
errorbar(muA,sigmaOff,'Color',reinfcolor+0.2,'LineWidth',2); hold on;
scatter([1,2,3],muA,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor); hold on;
errorbar(muB,sigmaOn,'Color',optocolor,'LineWidth',2);
scatter([1,2,3],muB,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
xticks([1,2,3]); xticklabels(conditions)
xlim([0.5 3.5]);box off;
ylabel('Hit Rate')
sigstar({[1,2]}, table2array(tbl2AnD(1,6)))
sigstar({[2,3]}, table2array(tbl2AnD(3,6)))
sigstar({[1,3]}, table2array(tbl2AnD(2,6)));ylim([0 1]);

%% subplot 3 Interaction plot for both groups FA rate

faRates=[allDataTestsOnly{6,39}(2,1:6)';allDataTestsOnly{6,40}(2,1:6)';allDataTestsOnly{6,41}(2,1:6)'];
optofaRates=[allDataTestsOnly{6,39}(4,1:6)';allDataTestsOnly{6,40}(4,1:6)';allDataTestsOnly{6,41}(4,1:6)'];
testFullTrial=horzcat(faRates, optofaRates);
sigmaOff    = [std(testFullTrial(1:6,1)), std(testFullTrial(7:12,1)), std(testFullTrial(13:18,1))];   % between-subject variability
sigmaOn    = [std(testFullTrial(1:6,2)), std(testFullTrial(7:12,2)), std(testFullTrial(13:18,2))];   % between-subject variability
% True means to create a main effect + interaction 
muA = [mean(testFullTrial(1:6,1)), mean(testFullTrial(7:12,1)), mean(testFullTrial(13:18,1))];    % Group A condition means
muB = [mean(testFullTrial(1:6,2)), mean(testFullTrial(7:12,2)), mean(testFullTrial(13:18,2))];    % Group B condition means (flatter -> interaction)
c2An = multcompare(anovaMat{5,3},"Estimate","row",'Display', 'off'); %full, stimulus, choicce
tbl2AnD = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
subplot(1,3,3);
errorbar(muA,sigmaOff,'Color',reinfcolor+0.2,'LineWidth',2); hold on;
scatter([1,2,3],muA,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor); hold on;
errorbar(muB,sigmaOn,'Color',optocolor,'LineWidth',2);
scatter([1,2,3],muB,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
xticks([1,2,3]); xticklabels(conditions);box off;
xlim([0.5 3.5]);ylim([0 1]);
ylabel('FA Rate');
sigstar({[1,2]}, table2array(tbl2AnD(1,6)))
sigstar({[2,3]}, table2array(tbl2AnD(3,6)))
sigstar({[1,3]}, table2array(tbl2AnD(2,6)))

interactionFig.Position(3:4)=[700 200];
saveas(gcf,['interactionPlotAnova']);
saveas(gcf,['interactionPlotAnova.png']);    
saveas(gcf,['interactionPlotAnova.pdf']);  

%% (Optional) Show the data in a table
% T = table(S, G, C, Y, 'VariableNames', {'Subject','Group','Condition','Y'});
% disp(T(1:12,:));