%To find GNG recordings that were made on the same day and bind the data of
%them in one output
%S.M. 11/27/2019
%S.M. revised 03/30/2020
function [dA]=bpod_filesbinder(matfiles2)

for ii=1:size(matfiles2,1) %from the matfiles given
    names(ii,:)=matfiles2(ii).name; %obtain the names of files
end
nA=str2num(names(:,end-18:end-11)); %only for LickGNG_Motiv experiment index the date of experiment (same day experiments will have == name)
if length(unique(nA))<size(matfiles2,1) %if there are repetitions of experiment
    [~,ind]=unique(nA); %obtain the indices of unique files
    duplicate_ind=setdiff(1:size(nA,1), ind); %obtain index of duplicated file
    duplicate_value=nA(duplicate_ind);
    udp=unique(duplicate_value);
    for iii=1:length(udp) %to get only the numbers that are duplicated in case more than one experiment was duplicated
        dA{iii}=find(nA==udp(iii)); %index of duplicated for each experiment
    end
else
    dA=NaN;
end
end


