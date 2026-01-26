function plotPsychCurve(allRpc,allDpc)
%% Psychometric curve (logistic) fit + plot
% Requires: Statistics and Machine Learning Toolbox (glmfit/glmval)

% vehicle
trials=50;
trialsfrac=trials*.01;
x   = [70 60 50 40]';     % stimulus levels (e.g., contrast)
k   = [mean(allRpc(:,1))*trialsfrac mean(allRpc(:,2))*trialsfrac mean(allRpc(:,3))*trialsfrac mean(allRpc(:,4))*trialsfrac]';              % # "correct" (or "yes") responses
n   = [trials trials trials trials]';              % total trials at each level

% Proportion correct/yes
p = k ./ n;

% --- Fit a logistic psychometric (binomial GLM with logit link) ---
% Model: logit(p) = b0 + b1*x
b = glmfit(x, [k n], 'binomial', 'link', 'logit');

% Dense x-grid for smooth curve
xq = linspace(min(x), max(x), 40)';
pq = glmval(b, xq, 'logit');  % predicted proportion

% --- (Optional) threshold at some criterion, e.g. 75% ---
crit = 0.75;
thr  = (log(crit/(1-crit)) - b(1)) / b(2);

figure;  hold on; box off
% Data points sized by trial count (optional)
scatter(x, p, 60, 'filled', 'MarkerFaceAlpha', 0.9); hold on;
% Fitted curve
plot(xq, pq, 'LineWidth', 2);

xlabel('Stimulus level');
ylabel('P(correct / yes)');
ylim([0.3 1]);
title('Psychometric curve');


% vehicle
x   = [70 60 50 40]';     % stimulus levels (e.g., contrast)
k   = [mean(allDpc(:,1))*trialsfrac mean(allDpc(:,2))*trialsfrac mean(allDpc(:,3))*trialsfrac mean(allDpc(:,4))*trialsfrac]';              % # "correct" (or "yes") responses
n   = [trials trials trials trials]';              % total trials at each level

% Proportion correct/yes
p = k ./ n;

% --- Fit a logistic psychometric (binomial GLM with logit link) ---
% Model: logit(p) = b0 + b1*x
b = glmfit(x, [k n], 'binomial', 'link', 'logit');

% Dense x-grid for smooth curve
xq = linspace(min(x), max(x), 400)';
pq = glmval(b, xq, 'logit');  % predicted proportion

% --- (Optional) threshold at some criterion, e.g. 75% ---
crit = 0.75;
thr  = (log(crit/(1-crit)) - b(1)) / b(2);

% Data points sized by trial count (optional)
scatter(x, p, 60, 'filled', 'MarkerFaceAlpha', 0.9);  hold on;

% Fitted curve
plot(xq, pq, 'LineWidth', 2);

xlabel('Stimulus level');
ylabel('P(correct / yes)');
ylim([0 1]);
title('Psychometric curve (logistic GLM)');
legend({'Vehicle Data', 'Vehicle Logistic fit','Ligand Data','Ligand Logistic fit'}, ...
       'Location','southeast');

% --- (Optional) 95% CI via parametric bootstrap ---
%{
B = 2000;
bh = mvnrnd(b, inv(b' * b), B); % quick-and-dirty; see note below
pq_boot = zeros(numel(xq), B);
for i = 1:B
    pq_boot(:,i) = glmval(bh(i,:)', xq, 'logit');
end
lo = prctile(pq_boot, 2.5, 2);
hi = prctile(pq_boot, 97.5, 2);
fill([xq; flipud(xq)], [lo; flipud(hi)], 1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
uistack(findobj(gca,'Type','line'),'top');
%}
