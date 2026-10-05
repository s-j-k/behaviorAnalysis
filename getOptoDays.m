function [days]=getOptoDays
    days={};
    days{1,1}='Animal';days{1,2}='ExpDays';days{1,3}='ExpRange';days{1,4}='Full Trial';days{1,5}='Stim';
    days{2,1}='sk163'; mgbDays = [0 0 0 0 1 1 0 0 1 0 1 0 0];
    mgbDays=logical(mgbDays);days{2,2}=mgbDays;expRange=1:11; days{2,3}=expRange;
    days{2,4} = [5,9]; days{2,5} = [6,11];
    
    days{3,1}='sk175'; mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0];mgbDays=logical(mgbDays);expRange=25:32;
    days{3,2}=mgbDays; days{3,3}=expRange;
    
    days{3,4}=[25,30];days{3,5}=[26,29];
    days{4,1}='sk176'; mgbDays=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 1];
    expRange=20:length(mgbDays);
    mgbDays=logical(mgbDays);days{4,2}=mgbDays; days{4,3}=expRange;
    days{4,4}=[23,25];%is this 35 or 25?
    days{4,5}=[26,31]; 
    days{5,1}='sk198'; mgbDays=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
        0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 1 0 0 0 0];
    mgbDays=logical(mgbDays);expRange=41:length(mgbDays); 
    days{5,2}=mgbDays; days{5,3}=expRange;
    days{5,4}=[41,44];    days{5,5}=[42,46];
    days{6,1}='sk203';mgbDays=[0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 0];
    mgbDays=logical(mgbDays); expRange=9:length(mgbDays); days{6,2}=mgbDays; days{6,3}=expRange;
    days{6,4}=[9,14];days{6,5}=[10,15];
    days{7,1}='sk204'; mgbDays=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 1 0];mgbDays=logical(mgbDays);
    expRange=27:length(mgbDays);days{7,2}=mgbDays; days{7,3}=expRange;
    days{7,4}=[27,32];days{7,5}=[30,34];   
    
    % need to populate this for the most recent cohort 13
    %these are animals 8 through 13 
    %cohort 14
    % iterate through getExpDays and get out the relevant information
    [subjlist, explist, pathsave]=getCohort(16);
    for ee=1:length(subjlist)
        [mgbDays, ~, expRange,~,subjNumber] = getExpDays(pathsave, ee);
        days{13+ee,1}=subjNumber;days{13+ee,2}=mgbDays;days{13+ee,2}=expRange;
        if length(expRange)>3
            days{13+ee,4}=[expRange(1),expRange(3)];days{13+ee,5}=[expRange(2),expRange(4)];
        elseif length(expRange)==2
            days{13+ee,4}=[expRange(1)];days{13+ee,5}=[expRange(2)];
        else
            days{13+ee,4}=[expRange(1)];
        end
    end
    
    
end
