function [pval]=mynormandttest(xwt,xko)

normthrwt=kstest(xwt);
normthrko=kstest(xko);
    if normthrwt==0 && normthrko==0
        [h pval]=ttest2(xwt,xko); %if data are normally distributed (0)
    else
        [pval h]=ranksum(xwt,xko); %if data are not normally distributed (1)
%         [pval h]=mwwtest(xwt,xko); %if data are not normally distributed (1)
    end
end 