function anovaInteractionPlot(aovMGBbyAn,statsMGBbyAn,aovIC,statsIC)
%% Repeated-measures example: 2 groups (between) x 3 conditions (within)
% Example data + "best practice" plots:
%  1) Spaghetti (each subject across conditions) within each group
%  2) Group mean ± within-subject 95% CI (Cousineau–Morey) overlaid

clear; close all; clc
rng(4);

groups     = ["Group A","Group B"];
conditions = ["C1","C2","C3"];

nSubPerGroup = 14;     % subjects per group
sigmaSubj    = 0.45;   % between-subject variability
sigmaWithin  = 0.25;   % within-subject noise

% True means to create a main effect + interaction (edit freely)
muA = [0.0, 0.3, 0.6];     % Group A condition means
muB = [0.1, 0.1, 0.1];     % Group B condition means (flatter -> interaction)

% --- Build long-format table-like arrays ---
Y = []; G = []; C = []; S = [];
subjCounter = 0;

for gi = 1:numel(groups)
    for si = 1:nSubPerGroup
        subjCounter = subjCounter + 1;
        subjID = subjCounter;

        subjOffset = sigmaSubj * randn; % subject random intercept

        if gi == 1
            mu = muA;
        else
            mu = muB;
        end

        y = mu(:) + subjOffset + sigmaWithin*randn(numel(conditions),1);

        Y = [Y; y];
        G = [G; repmat(groups(gi), numel(conditions), 1)];
        C = [C; conditions(:)];
        S = [S; repmat(subjID, numel(conditions), 1)];
    end
end

%% --- Helper: compute group means + within-subject 95% CI (Cousineau–Morey) ---
% Steps:
% 1) Normalize within each subject (remove subject mean, add grand mean)
% 2) Compute SEM of normalized values per condition within each group
% 3) Apply Morey correction: sqrt(k/(k-1)) where k=#conditions
k = numel(conditions);
morey = sqrt(k/(k-1));

m  = nan(numel(groups), k);
lo = nan(size(m));
hi = nan(size(m));

for gi = 1:numel(groups)
    % Grab this group's data
    idxG = (G == groups(gi));
    Yg = Y(idxG);
    Cg = C(idxG);
    Sg = S(idxG);

    subjIDs = unique(Sg);
    Ynorm = nan(size(Yg));

    grandMean = mean(Yg);

    % Normalize per subject
    for s = subjIDs'
        idxS = (Sg == s);
        subjMean = mean(Yg(idxS));
        Ynorm(idxS) = Yg(idxS) - subjMean + grandMean;
    end

    % Compute mean and CI per condition
    for ci = 1:k
        idxC = (Cg == conditions(ci));
        yRaw  = Yg(idxC);
        yNorm = Ynorm(idxC);

        m(gi,ci) = mean(yRaw);

        n = numel(subjIDs); % subjects (not trials)
        sem = std(yNorm,0) / sqrt(n);
        sem = sem * morey;

        tcrit = tinv(0.975, n-1);
        lo(gi,ci) = m(gi,ci) - tcrit*sem;
        hi(gi,ci) = m(gi,ci) + tcrit*sem;
    end
end

%% =========================
% Figure 1: Spaghetti per group (paired lines) + mean ± within-subject CI
% =========================
figure('Color','w','Name','Repeated measures: spaghetti + mean/CI');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

x = 1:k;

for gi = 1:numel(groups)
    nexttile; hold on; box on

    % Plot each subject line (spaghetti)
    idxG = (G == groups(gi));
    subjIDs = unique(S(idxG));

    for s = subjIDs'
        ys = nan(1,k);
        for ci = 1:k
            idx = idxG & (S==s) & (C==conditions(ci));
            ys(ci) = Y(idx);
        end
        plot(x, ys, '-', 'LineWidth', 0.9); % many lines; keep thin
        scatter(x, ys, 22, 'filled', 'MarkerFaceAlpha', 0.85);
    end

    % Overlay mean ± within-subject CI
    errorbar(x, m(gi,:), m(gi,:)-lo(gi,:), hi(gi,:)-m(gi,:), ...
             'k-o', 'LineWidth', 2.2, 'CapSize', 10, 'MarkerSize', 7);

    xlim([0.7, k+0.3])
    xticks(x); xticklabels(conditions)
    ylabel('Dependent variable')
    title(sprintf('%s', groups(gi)))
end

sgtitle('Repeated-measures visualization: subject trajectories + mean \pm within-subject 95% CI')

%% =========================
% Figure 2: Interaction plot (means ± within-subject CI) for both groups
% =========================
figure('Color','w','Name','Repeated measures: interaction (mean/CI)'); hold on; box on
for gi = 1:numel(groups)
    errorbar(x, m(gi,:), m(gi,:)-lo(gi,:), hi(gi,:)-m(gi,:), ...
             '-o', 'LineWidth', 2.2, 'CapSize', 10, 'MarkerSize', 7);
end
xticks(x); xticklabels(conditions)
ylabel('Mean dependent variable')
title('Interaction plot (mean \pm within-subject 95% CI)')
legend(groups, 'Location','best')

%% (Optional) Show the example data in a table
% T = table(S, G, C, Y, 'VariableNames', {'Subject','Group','Condition','Y'});
% disp(T(1:12,:));