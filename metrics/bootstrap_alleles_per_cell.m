function [apc_mu, apc_sig] = bootstrap_alleles_per_cell(summary, N_trials)
    assert(isa(summary, 'ExperimentSummary'));
    N_alleles = size(summary.alleles,1);
    N_observations = sum(summary.allele_freqs);
    apc_mu = N_alleles/N_observations;
    datavec = repelem([1:N_alleles]', summary.allele_freqs);
    apc_trial = zeros(N_trials,1);
    for i = 1:N_trials
        apc_trial(i) = length(unique(datasample(datavec, N_observations, 'Replace', true)))/N_observations;
    end
    apc_sig = sqrt(sum((apc_trial-apc_mu).^2)/N_trials);
end
