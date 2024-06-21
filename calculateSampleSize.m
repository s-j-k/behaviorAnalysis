function sample_size = calculateSampleSize(effect_size, alpha, power, threshold)
    if nargin < 4; threshold = 0.001; end
    effect_size = abs(effect_size); % to avoid inf ss_estimate
    ss_estimate = 30; % initial sample size estimate
    while true
        df                   = 2*(ss_estimate - 1);
        critical_value       = tinv(1 - alpha/2, df);
        non_centrality_param = effect_size*sqrt(ss_estimate/2);
        current_power        = 1 - nctcdf(critical_value, df, non_centrality_param);

        if abs(current_power - power) < threshold
            break;
        end
        
        ss_estimate = ss_estimate*(power/current_power);
    end
    
    sample_size = ss_estimate;
end