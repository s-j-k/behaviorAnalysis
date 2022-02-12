function [rho,yfit] = LinearRegression(x,y)

p = polyfit(x,y,1);% gives slope and intercept of linear predictor 
yfit = polyval(p,x); % use polyval to predict y given the polyfit estimate 'p' 
yresid = y - yfit; % Compute the residual values as a vector signed numbers: 
SSresid = sum(yresid.^2); % Square the residuals and total them obtain the residual sum of squares: 
SStotal = (length(y)-1)*var(y); % Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1: 
rho = sign(p(1))*(1-SSresid/SStotal); % Compute R²
