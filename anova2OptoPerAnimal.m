function [aovMGB,statsMGB,aovIC,statsIC]=anova2OptoPerAnimal(mgbTempTestsOnly,tempTestsOnly)

% now per animal statistical analysis as a 2way anova
% no light vs. light all inactivation condition
% MGB
testFullTrial=NaN(18,2); 
% Full Trial, light OFF column 1, light ON column 2
testFullTrial(1:2,1)=mgbTempTestsOnly{2,27}(~isnan(mgbTempTestsOnly{2,27}));
testFullTrial(1:2,2)=mgbTempTestsOnly{2,28}(~isnan(mgbTempTestsOnly{2,28}));
testFullTrial(3:4,1)=mgbTempTestsOnly{3,27}(~isnan(mgbTempTestsOnly{3,27}));
testFullTrial(3:4,2)=mgbTempTestsOnly{3,28}(~isnan(mgbTempTestsOnly{3,28}));
fullTrialtemp=mgbTempTestsOnly{4,27}(~isnan(mgbTempTestsOnly{4,27}));
testFullTrial(5:6,1)=fullTrialtemp(1:2);
fullTrialtemp=mgbTempTestsOnly{4,28}(~isnan(mgbTempTestsOnly{4,28}));
testFullTrial(5:6,2)=fullTrialtemp(1:2);
% tone Trial, light OFF column 1, light ON column 2
testFullTrial(7:8,1)=mgbTempTestsOnly{2,29}(~isnan(mgbTempTestsOnly{2,29}));
testFullTrial(7:8,2)=mgbTempTestsOnly{2,30}(~isnan(mgbTempTestsOnly{2,30}));
testFullTrial(9:10,1)=mgbTempTestsOnly{3,29}(~isnan(mgbTempTestsOnly{3,29}));
testFullTrial(9:10,2)=mgbTempTestsOnly{3,30}(~isnan(mgbTempTestsOnly{3,30}));
fullTrialtemp=mgbTempTestsOnly{4,29}(~isnan(mgbTempTestsOnly{4,29}));
testFullTrial(11:12,1)=fullTrialtemp(1:2);
fullTrialtemp=mgbTempTestsOnly{4,30}(~isnan(mgbTempTestsOnly{4,30}));
testFullTrial(11:12,2)=fullTrialtemp(1:2);
% Choice Trial, light OFF column 1, light ON column 2
fullTrialtemp=mgbTempTestsOnly{2,31}(~isnan(mgbTempTestsOnly{2,31}));
testFullTrial(13:14,1)=fullTrialtemp(1:2);
fullTrialtemp=mgbTempTestsOnly{2,32}(~isnan(mgbTempTestsOnly{2,32}));
testFullTrial(13:14,2)=fullTrialtemp(1:2);
fullTrialtemp=mgbTempTestsOnly{3,31}(~isnan(mgbTempTestsOnly{3,31}));
testFullTrial(15:16,1)=fullTrialtemp(1:2);
fullTrialtemp=mgbTempTestsOnly{3,32}(~isnan(mgbTempTestsOnly{3,32}));
testFullTrial(15:16,2)=fullTrialtemp(1:2);
fullTrialtemp=mgbTempTestsOnly{4,31}(~isnan(mgbTempTestsOnly{4,31}));
testFullTrial(17:18,1)=fullTrialtemp(1:2);
fullTrialtemp=mgbTempTestsOnly{4,32}(~isnan(mgbTempTestsOnly{4,32}));
testFullTrial(17:18,2)=fullTrialtemp(1:2);

[aovMGB,~,statsMGB]=anova2(testFullTrial,6); % For the MGB 
c1 = multcompare(statsMGB);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2 = multcompare(statsMGB,"Estimate","row");
tbl2 = array2table(c2,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% figure; % i was trying to plot something here? not sure. 
% light = [repmat("On",18,1);repmat("Off",18,1);];
% optoCond = [repmat("Full Trial",6,1);repmat("Tone",6,1);repmat("Choice",6,1)];
% factors = {light, optoCond};
% aov = anova(factors,testFullTrial(:),FactorNames=["Light","optoCond"],...
%     ModelSpecification="interactions");
% boxchart(aov,["Light on vs. Off","Opto Condition"])


% now do the IC
testFullTrial=NaN(22,2); 
% Full Trial, light OFF column 1, light ON column 2
testFullTrial(1:4,1)=tempTestsOnly{2,33}(~isnan(tempTestsOnly{2,33}));
testFullTrial(1:4,2)=tempTestsOnly{2,34}(~isnan(tempTestsOnly{2,34}));
testFullTrial(5:8,1)=tempTestsOnly{3,33}(~isnan(tempTestsOnly{3,33}));
testFullTrial(5:8,2)=tempTestsOnly{3,34}(~isnan(tempTestsOnly{3,34}));
testFullTrial(1,:)=[]; % to have 7 samples for each
% tone Trial, light OFF column 1, light ON column 2
testFullTrial(8:10,1)=tempTestsOnly{2,35}(~isnan(tempTestsOnly{2,35}));
testFullTrial(8:10,2)=tempTestsOnly{2,36}(~isnan(tempTestsOnly{2,36}));
testFullTrial(11:14,1)=tempTestsOnly{3,35}(~isnan(tempTestsOnly{3,35}));
testFullTrial(11:14,2)=tempTestsOnly{3,36}(~isnan(tempTestsOnly{3,36}));
% Choice Trial, light OFF column 1, light ON column 2
fullTrialtemp=tempTestsOnly{2,37}(~isnan(tempTestsOnly{2,37}));
testFullTrial(15:17,1)=fullTrialtemp(1:3);
fullTrialtemp=tempTestsOnly{2,38}(~isnan(tempTestsOnly{2,38}));
testFullTrial(15:17,2)=fullTrialtemp(1:3);
fullTrialtemp=tempTestsOnly{3,37}(~isnan(tempTestsOnly{3,37}));
testFullTrial(18:21,1)=fullTrialtemp(1:4);
fullTrialtemp=tempTestsOnly{3,38}(~isnan(tempTestsOnly{3,38}));
testFullTrial(18:21,2)=fullTrialtemp(1:4);

[aovIC,~,statsIC]=anova2(testFullTrial,7); 
c1IC = multcompare(statsIC);
tbl1IC = array2table(c1IC,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2IC = multcompare(statsIC,"Estimate","row");
tbl2 = array2table(c2IC,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);