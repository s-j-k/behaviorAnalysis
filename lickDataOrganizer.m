function [lickMat]=lickDataOrganizer(allDataTestsOnly)
%% organize lick data
lickMat={};
optoflag=1;subjlist=NaN;
for nbsubj=2:size(allDataTestsOnly,1)
    lickDays=days{nbsubj,2};
    expRange=days{nbsubj,3};lickDays=lickDays(expRange);
    matrix=allDataTestsOnly{nbsubj,26};
    lickMatOut=getLickLatHist(matrix,nbsubj,subjlist,lickDays,expRange,optoflag);
    lickMat{nbsubj}=lickMatOut;
end
end