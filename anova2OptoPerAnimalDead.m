function [aovMGB,statsMGB]=anova2OptoPerAnimalDead(allDataTestsOnly)

% now per animal statistical analysis as a 2way anova
% no light vs. light all inactivation condition
% MGB
testFullTrial=NaN(30,2); 
% Dead 1
testFullTrial(1:2,1)=allDataTestsOnly{2,19};
testFullTrial(1:2,2)=allDataTestsOnly{2,20};
testFullTrial(3:4,1)=allDataTestsOnly{3,19};
testFullTrial(3:4,2)=allDataTestsOnly{3,20};
testFullTrial(5:6,1)=allDataTestsOnly{4,19};
testFullTrial(5:6,2)=allDataTestsOnly{4,20}; 
% Dead 2
testFullTrial(7:8,1)=allDataTestsOnly{2,21};
testFullTrial(7:8,2)=allDataTestsOnly{2,22};
testFullTrial(9:10,1)=allDataTestsOnly{3,21};
testFullTrial(9:10,2)=allDataTestsOnly{3,22};
testFullTrial(11:12,1)=allDataTestsOnly{4,21};
testFullTrial(11:12,2)=allDataTestsOnly{4,22}; 
% Dead 3
testFullTrial(13:14,1)=allDataTestsOnly{2,23};
testFullTrial(13:14,2)=allDataTestsOnly{2,24};
testFullTrial(15:16,1)=allDataTestsOnly{3,23};
testFullTrial(15:16,2)=allDataTestsOnly{3,24};
testFullTrial(17:18,1)=allDataTestsOnly{4,23};
testFullTrial(17:18,2)=allDataTestsOnly{4,24}; 
% Dead 4
testFullTrial(19:20,1)=allDataTestsOnly{2,25};
testFullTrial(19:20,2)=allDataTestsOnly{2,26};
testFullTrial(21:22,1)=allDataTestsOnly{3,25};
testFullTrial(21:22,2)=allDataTestsOnly{3,26};
testFullTrial(23:24,1)=allDataTestsOnly{4,25};
testFullTrial(23:24,2)=allDataTestsOnly{4,26}; 
% Dead 5
testFullTrial(25:26,1)=allDataTestsOnly{2,27};
testFullTrial(25:26,2)=allDataTestsOnly{2,28};
testFullTrial(27:28,1)=allDataTestsOnly{3,27};
testFullTrial(27:28,2)=allDataTestsOnly{3,28};
testFullTrial(29:30,1)=allDataTestsOnly{4,27};
testFullTrial(29:30,2)=allDataTestsOnly{4,28}; 
[aovMGB,~,statsMGB]=anova2(testFullTrial,6); % For the MGB 
c1 = multcompare(statsMGB);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2 = multcompare(statsMGB,"Estimate","row");
tbl2 = array2table(c2,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
