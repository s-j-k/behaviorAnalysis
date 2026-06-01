function [anovaMat]=anova2OptoPerAnimal(allDataTestsOnly)

% per session 2way anova
% no light vs. light all inactivation condition
% MGB by session (two sessions per condition, per animal
testFullTrial=NaN(36,2); 
% % Full Trial, light OFF column 1, light ON column 2
testFullTrial(1:2,1)=allDataTestsOnly{2,27}(~isnan(allDataTestsOnly{2,27}));
testFullTrial(1:2,2)=allDataTestsOnly{2,28}(~isnan(allDataTestsOnly{2,28}));
testFullTrial(3:4,1)=allDataTestsOnly{3,27}(~isnan(allDataTestsOnly{3,27}));
testFullTrial(3:4,2)=allDataTestsOnly{3,28}(~isnan(allDataTestsOnly{3,28}));
fullTrialtemp=allDataTestsOnly{4,27}(~isnan(allDataTestsOnly{4,27}));
testFullTrial(5:6,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,28}(~isnan(allDataTestsOnly{4,28}));
testFullTrial(5:6,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{5,27}(~isnan(allDataTestsOnly{5,27}));
testFullTrial(7:8,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{5,28}(~isnan(allDataTestsOnly{5,28}));
testFullTrial(7:8,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{6,27}(~isnan(allDataTestsOnly{6,27}));
testFullTrial(9:10,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{6,28}(~isnan(allDataTestsOnly{6,28}));
testFullTrial(9:10,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{7,27}(~isnan(allDataTestsOnly{7,27}));
testFullTrial(11:12,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{7,28}(~isnan(allDataTestsOnly{7,28}));
testFullTrial(11:12,2)=fullTrialtemp(1:2);
% tone Trial, light OFF column 1, light ON column 2
testFullTrial(13:14,1)=allDataTestsOnly{2,29}(~isnan(allDataTestsOnly{2,29}));
testFullTrial(13:14,2)=allDataTestsOnly{2,30}(~isnan(allDataTestsOnly{2,30}));
testFullTrial(15:16,1)=allDataTestsOnly{3,29}(~isnan(allDataTestsOnly{3,29}));
testFullTrial(15:16,2)=allDataTestsOnly{3,30}(~isnan(allDataTestsOnly{3,30}));
fullTrialtemp=allDataTestsOnly{4,29}(~isnan(allDataTestsOnly{4,29}));
testFullTrial(17:18,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,30}(~isnan(allDataTestsOnly{4,30}));
testFullTrial(17:18,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{5,29}(~isnan(allDataTestsOnly{5,29}));
testFullTrial(19:20,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{5,30}(~isnan(allDataTestsOnly{5,30}));
testFullTrial(19:20,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{6,29}(~isnan(allDataTestsOnly{6,29}));
testFullTrial(21:22,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{6,30}(~isnan(allDataTestsOnly{6,30}));
testFullTrial(21:22,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{7,29}(~isnan(allDataTestsOnly{7,29}));
testFullTrial(23:24,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{7,30}(~isnan(allDataTestsOnly{7,30}));
testFullTrial(23:24,2)=fullTrialtemp(1:2);
% Choice Trial, light OFF column 1, light ON column 2
fullTrialtemp=allDataTestsOnly{2,31}(~isnan(allDataTestsOnly{2,31}));
testFullTrial(25:26,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{2,32}(~isnan(allDataTestsOnly{2,32}));
testFullTrial(25:26,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{3,31}(~isnan(allDataTestsOnly{3,31}));
testFullTrial(27:28,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{3,32}(~isnan(allDataTestsOnly{3,32}));
testFullTrial(27:28,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,31}(~isnan(allDataTestsOnly{4,31}));
testFullTrial(29:30,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{4,32}(~isnan(allDataTestsOnly{4,32}));
testFullTrial(29:30,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{5,31}(~isnan(allDataTestsOnly{5,31}));
testFullTrial(31:32,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{5,32}(~isnan(allDataTestsOnly{5,32}));
testFullTrial(31:32,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{6,31}(~isnan(allDataTestsOnly{6,31}));
testFullTrial(33:34,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{6,32}(~isnan(allDataTestsOnly{6,32}));
testFullTrial(33:34,2)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{7,31}(~isnan(allDataTestsOnly{7,31}));
testFullTrial(35:36,1)=fullTrialtemp(1:2);
fullTrialtemp=allDataTestsOnly{7,32}(~isnan(allDataTestsOnly{7,32}));
testFullTrial(35:36,2)=fullTrialtemp(1:2);
[aovMGBbySess,~,statsMGBbySess]=anova2(testFullTrial,12); % For the MGB 
c1 = multcompare(statsMGBbySess);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2 = multcompare(statsMGBbySess,"Estimate","row");
tbl2 = array2table(c2,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

% now per animal
testFullTrialByAnimal=NaN(18,2); 
% Full Trial, light OFF column 1, light ON column 2
testFullTrialByAnimal(1,1)=mean(allDataTestsOnly{2,27}(~isnan(allDataTestsOnly{2,27})));
testFullTrialByAnimal(1,2)=mean(allDataTestsOnly{2,28}(~isnan(allDataTestsOnly{2,28})));
testFullTrialByAnimal(2,1)=mean(allDataTestsOnly{3,27}(~isnan(allDataTestsOnly{3,27})));
testFullTrialByAnimal(2,2)=mean(allDataTestsOnly{3,28}(~isnan(allDataTestsOnly{3,28})));
fullTrialtemp=allDataTestsOnly{4,27}(~isnan(allDataTestsOnly{4,27}));
testFullTrialByAnimal(3,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{4,28}(~isnan(allDataTestsOnly{4,28}));
testFullTrialByAnimal(3,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{5,27}(~isnan(allDataTestsOnly{5,27}));
testFullTrialByAnimal(4,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{5,28}(~isnan(allDataTestsOnly{5,28}));
testFullTrialByAnimal(4,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{6,27}(~isnan(allDataTestsOnly{6,27}));
testFullTrialByAnimal(5,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{6,28}(~isnan(allDataTestsOnly{6,28}));
testFullTrialByAnimal(5,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{7,27}(~isnan(allDataTestsOnly{7,27}));
testFullTrialByAnimal(6,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{7,28}(~isnan(allDataTestsOnly{7,28}));
testFullTrialByAnimal(6,2)=mean(fullTrialtemp(1:2));
% tone Trial, light OFF column 1, light ON column 2
testFullTrialByAnimal(7,1)=mean(allDataTestsOnly{2,29}(~isnan(allDataTestsOnly{2,29})));
testFullTrialByAnimal(7,2)=mean(allDataTestsOnly{2,30}(~isnan(allDataTestsOnly{2,30})));
testFullTrialByAnimal(8,1)=mean(allDataTestsOnly{3,29}(~isnan(allDataTestsOnly{3,29})));
testFullTrialByAnimal(8,2)=mean(allDataTestsOnly{3,30}(~isnan(allDataTestsOnly{3,30})));
fullTrialtemp=allDataTestsOnly{4,29}(~isnan(allDataTestsOnly{4,29}));
testFullTrialByAnimal(9,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{4,30}(~isnan(allDataTestsOnly{4,30}));
testFullTrialByAnimal(9,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{5,29}(~isnan(allDataTestsOnly{5,29}));
testFullTrialByAnimal(10,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{5,30}(~isnan(allDataTestsOnly{5,30}));
testFullTrialByAnimal(10,2)=mean(fullTrialtemp(1:2));

fullTrialtemp=allDataTestsOnly{6,29}(~isnan(allDataTestsOnly{6,29}));
testFullTrialByAnimal(11,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{6,30}(~isnan(allDataTestsOnly{6,30}));
testFullTrialByAnimal(11,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{7,29}(~isnan(allDataTestsOnly{7,29}));
testFullTrialByAnimal(12,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{7,30}(~isnan(allDataTestsOnly{7,30}));
testFullTrialByAnimal(12,2)=mean(fullTrialtemp(1:2));
% Choice Trial, light OFF column 1, light ON column 2
fullTrialtemp=allDataTestsOnly{2,31}(~isnan(allDataTestsOnly{2,31}));
testFullTrialByAnimal(13,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{2,32}(~isnan(allDataTestsOnly{2,32}));
testFullTrialByAnimal(13,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{3,31}(~isnan(allDataTestsOnly{3,31}));
testFullTrialByAnimal(14,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{3,32}(~isnan(allDataTestsOnly{3,32}));
testFullTrialByAnimal(14,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{4,31}(~isnan(allDataTestsOnly{4,31}));
testFullTrialByAnimal(15,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{4,32}(~isnan(allDataTestsOnly{4,32}));
testFullTrialByAnimal(15,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{5,31}(~isnan(allDataTestsOnly{5,31}));
testFullTrialByAnimal(16,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{5,32}(~isnan(allDataTestsOnly{5,32}));
testFullTrialByAnimal(16,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{6,31}(~isnan(allDataTestsOnly{6,31}));
testFullTrialByAnimal(17,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{6,32}(~isnan(allDataTestsOnly{6,32}));
testFullTrialByAnimal(17,2)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{7,31}(~isnan(allDataTestsOnly{7,31}));
testFullTrialByAnimal(18,1)=mean(fullTrialtemp(1:2));
fullTrialtemp=allDataTestsOnly{7,32}(~isnan(allDataTestsOnly{7,32}));
testFullTrialByAnimal(18,2)=mean(fullTrialtemp(1:2));
[aovMGBbyAn,~,statsMGBbyAn]=anova2(testFullTrialByAnimal,6); % For the MGB 
c1An = multcompare(statsMGBbyAn);
tbl1An = array2table(c1An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2An = multcompare(statsMGBbyAn,"Estimate","row");
tbl2An = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

% % now do the IC
% testFullTrial=NaN(22,2); 
% % Full Trial, light OFF column 1, light ON column 2
% testFullTrial(1:4,1)=tempTestsOnly{2,33}(~isnan(tempTestsOnly{2,33}));
% testFullTrial(1:4,2)=tempTestsOnly{2,34}(~isnan(tempTestsOnly{2,34}));
% testFullTrial(5:8,1)=tempTestsOnly{3,33}(~isnan(tempTestsOnly{3,33}));
% testFullTrial(5:8,2)=tempTestsOnly{3,34}(~isnan(tempTestsOnly{3,34}));
% testFullTrial(1,:)=[]; % to have 7 samples for each
% % tone Trial, light OFF column 1, light ON column 2
% testFullTrial(8:10,1)=tempTestsOnly{2,35}(~isnan(tempTestsOnly{2,35}));
% testFullTrial(8:10,2)=tempTestsOnly{2,36}(~isnan(tempTestsOnly{2,36}));
% testFullTrial(11:14,1)=tempTestsOnly{3,35}(~isnan(tempTestsOnly{3,35}));
% testFullTrial(11:14,2)=tempTestsOnly{3,36}(~isnan(tempTestsOnly{3,36}));
% % Choice Trial, light OFF column 1, light ON column 2
% fullTrialtemp=tempTestsOnly{2,37}(~isnan(tempTestsOnly{2,37}));
% testFullTrial(15:17,1)=fullTrialtemp(1:3);
% fullTrialtemp=tempTestsOnly{2,38}(~isnan(tempTestsOnly{2,38}));
% testFullTrial(15:17,2)=fullTrialtemp(1:3);
% fullTrialtemp=tempTestsOnly{3,37}(~isnan(tempTestsOnly{3,37}));
% testFullTrial(18:21,1)=fullTrialtemp(1:4);
% fullTrialtemp=tempTestsOnly{3,38}(~isnan(tempTestsOnly{3,38}));
% testFullTrial(18:21,2)=fullTrialtemp(1:4);
% 
% [aovIC,~,statsIC]=anova2(testFullTrial,7); 
% c1IC = multcompare(statsIC);
% tbl1IC = array2table(c1IC,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% c2IC = multcompare(statsIC,"Estimate","row");
% tbl2 = array2table(c2IC,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

%% d' stats
testFullTrial=NaN(18,2); 
testFullTrial(1:6,1:2)=allDataTestsOnly{2,39}';
testFullTrial(7:12,1:2)=allDataTestsOnly{2,40}';
testFullTrial(13:18,1:2)=allDataTestsOnly{2,41}';
[aovMGBbyAnDprime,~,statsMGBbyAnDprime]=anova2(testFullTrial,6);
c1An = multcompare(statsMGBbyAnDprime);
tbl1An = array2table(c1An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2An = multcompare(statsMGBbyAnDprime,"Estimate","row");
tbl2AnD = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
mean(allDataTestsOnly{2,40}(1,:));
std(allDataTestsOnly{2,40}(1,:));
mean(allDataTestsOnly{2,40}(2,:));
std(allDataTestsOnly{2,40}(2,:));

hitStats=NaN(18,2); 
hitStats(1:6,1)=allDataTestsOnly{6,39}(1,1:6)';
hitStats(1:6,2)=allDataTestsOnly{6,39}(3,1:6)';
hitStats(7:12,1)=allDataTestsOnly{6,40}(1,1:6)';
hitStats(7:12,2)=allDataTestsOnly{6,40}(3,1:6)';
hitStats(13:18,1)=allDataTestsOnly{6,41}(1,1:6)';
hitStats(13:18,2)=allDataTestsOnly{6,41}(3,1:6)';
[aovMGBbyAnHit,~,statsMGBbyAnHit]=anova2(hitStats,6);
c1An = multcompare(statsMGBbyAnHit);
tbl1An = array2table(c1An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2An = multcompare(statsMGBbyAnHit,"Estimate","row");
tbl2AnHit = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

FAStats=NaN(18,2); 
FAStats(1:6,1)=allDataTestsOnly{6,39}(2,1:6)';
FAStats(1:6,2)=allDataTestsOnly{6,39}(4,1:6)';
FAStats(7:12,1)=allDataTestsOnly{6,40}(2,1:6)';
FAStats(7:12,2)=allDataTestsOnly{6,40}(4,1:6)';
FAStats(13:18,1)=allDataTestsOnly{6,41}(2,1:6)';
FAStats(13:18,2)=allDataTestsOnly{6,41}(4,1:6)';
[aovMGBbyAnFA,~,statsMGBbyAnFA]=anova2(FAStats,6);
c1An = multcompare(statsMGBbyAnFA);
tbl1An = array2table(c1An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2An = multcompare(statsMGBbyAnFA,"Estimate","row");
tbl2AnFA = array2table(c2An,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);


anovaMat={};
anovaMat{1,2}='p values';anovaMat{1,3}='stat values';
anovaMat{2,1}='By Session percent correct';anovaMat{2,2}=aovMGBbySess;anovaMat{2,3}=statsMGBbySess;
anovaMat{3,1}='By animal percent correct';anovaMat{3,2}=aovMGBbyAn;anovaMat{3,3}=statsMGBbyAn;
anovaMat{4,1}='By Animal Hit Rate';anovaMat{4,2}=aovMGBbyAnHit;anovaMat{4,3}=statsMGBbyAnHit;
anovaMat{5,1}='By Animal FA Rate';anovaMat{5,2}=aovMGBbyAnFA;anovaMat{5,3}=statsMGBbyAnFA;
anovaMat{6,1}='By Animal MGB Dprime';anovaMat{6,2}=aovMGBbyAnDprime;anovaMat{6,3}=statsMGBbyAnDprime;
% anovaMat{7,1}='By Animal IC Percent Correct';anovaMat{7,2}=aovIC;anovaMat{7,3}=statsIC;
end
