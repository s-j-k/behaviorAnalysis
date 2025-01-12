function [aovMGB,statsMGB]=anova2OptoPerAnimalDead(allDataTestsOnly)

% now per animal statistical analysis as a 2way anova
% no light vs. light all inactivation condition
% MGB
testFullTrial=NaN(25,2); 
% Dead 1
testFullTrial(1:2,1)=allDataTestsOnly{2,19};
testFullTrial(1:2,2)=allDataTestsOnly{2,20};
testFullTrial(3,1)=allDataTestsOnly{3,19}(1,1);
testFullTrial(3,2)=allDataTestsOnly{3,20}(1,1);
testFullTrial(4:5,1)=allDataTestsOnly{4,19};
testFullTrial(4:5,2)=allDataTestsOnly{4,20}; 
% Dead 2
testFullTrial(6:7,1)=allDataTestsOnly{2,21};
testFullTrial(6:7,2)=allDataTestsOnly{2,22};
testFullTrial(8,1)=allDataTestsOnly{3,21}(1,1);
testFullTrial(8,2)=allDataTestsOnly{3,22}(1,1);
testFullTrial(9:10,1)=allDataTestsOnly{4,21};
testFullTrial(9:10,2)=allDataTestsOnly{4,22}; 
% Dead 3
testFullTrial(11:12,1)=allDataTestsOnly{2,23};
testFullTrial(11:12,2)=allDataTestsOnly{2,24};
testFullTrial(13,1)=allDataTestsOnly{3,23}(1,1);
testFullTrial(13,2)=allDataTestsOnly{3,24}(1,1);
testFullTrial(14:15,1)=allDataTestsOnly{4,23};
testFullTrial(14:15,2)=allDataTestsOnly{4,24}; 
% Dead 4
testFullTrial(16:17,1)=allDataTestsOnly{2,25};
testFullTrial(16:17,2)=allDataTestsOnly{2,26};
testFullTrial(18,1)=allDataTestsOnly{3,25}(1,1);
testFullTrial(18,2)=allDataTestsOnly{3,26}(1,1);
testFullTrial(19:20,1)=allDataTestsOnly{4,25};
testFullTrial(19:20,2)=allDataTestsOnly{4,26}; 
% Dead 5
testFullTrial(21:22,1)=allDataTestsOnly{2,27};
testFullTrial(21:22,2)=allDataTestsOnly{2,28};
testFullTrial(23,1)=allDataTestsOnly{3,27}(1,1);
testFullTrial(23,2)=allDataTestsOnly{3,28}(1,1);
testFullTrial(24:25,1)=allDataTestsOnly{4,27};
testFullTrial(24:25,2)=allDataTestsOnly{4,28}; 
[aovMGB,~,statsMGB]=anova2(testFullTrial,5); % For the MGB 
c1 = multcompare(statsMGB);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
c2 = multcompare(statsMGB,"Estimate","row");
tbl2 = array2table(c2,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
