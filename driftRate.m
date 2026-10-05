function v = driftRate(t, ~, p, evidence, laserCondition, toneDuration)

    tonePeriod   = (t >= 0) && (t < toneDuration);
    choicePeriod = (t >= toneDuration);

    % Baseline drift, including retained evidence after tone offset
    v = p.beta0 ...
      + p.betaEvidenceTone     * evidence * tonePeriod ...
      + p.betaRetainedEvidence * evidence * choicePeriod;

    switch laserCondition
        case 0
            % Laser off
            laserActive = false;
            conditionIndex = [];

        case 1
            % Full-trial inactivation
            laserActive = true;
            conditionIndex = 1;

        case 2
            % stimulus-period inactivation
            laserActive = tonePeriod;
            conditionIndex = 2;

        case 3
            % Choice-period inactivation
            laserActive = choicePeriod;
            conditionIndex = 3;

        otherwise
            error("Unknown laser condition: %g", laserCondition);
    end

    if laserActive
        v = v ...
          + p.betaLaser(conditionIndex) ...
          + p.betaLaserEvidence(conditionIndex) * evidence;
    end
end