function [fitresult, gof]=mysigmfit(xf,yf)

clear fitresult gof opts ft af bf cf df
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.MaxIter=1000;
opts.MaxFunEvals=1000;
ft=fittype('a+(b-a)/(1+exp(-(x-c)/d))');
% a is the lower bound, starting d' 
% b is the upper bound, ending d'
% c is halfway between a and b, in units of the x axis values--in a usual
% sigmoid this is the 0.5 point (since sigmoids are usually standardized to
% 0 to 1
% d is the slope of the sigmoid, in units of the x axis--inverse slope 
% another way to think of d is the number of trials the animal takes to
% learn 

%in case there are NaNs in the vectors
xf(isnan(yf))=[];
yf(isnan(yf))=[]; 

af=mean(yf(1:3)); %min value
bf=mean(yf(end-2:end)); %max value
[vif iif]=min(abs(yf-[(bf-af)/2])); %find middle point in y and index for x
cf=xf(iif); %inflection point
df=(max(xf)-min(xf))/(bf-af); %slope
opts.StartPoint=[af bf cf df];

[fitresult, gof] = fit(xf, yf, ft, opts);
end

    
%prev values
%     bf=mean(yf(1:3)); %min value
%     cf=mean(yf(end-2:end)); %max value
%     af=xf(iif); %inflection point
%     df=1; %slope
%     [vif iif]=min(abs(yf-[(cf-bf)/2])); %find middle point in y and index for x
%     ft =  fittype('c + (b - c)/(1 + (x/a)^d)'); %fit function
