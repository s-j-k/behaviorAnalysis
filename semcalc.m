function sems=semcalc(y)

num=sum(~isnan(y));
sems=nanstd(y)./sqrt(num);
end 