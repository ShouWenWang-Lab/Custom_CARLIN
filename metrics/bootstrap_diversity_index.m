function [di_mu, di_sig, di_trial] = bootstrap_diversity_index(summary, N_trials)
    assert(isa(summary, 'ExperimentSummary'));
    refseq = CARLIN_def.getInstance.seq.CARLIN;
    is_template = ismember(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), refseq);
    N_total  = sum(summary.allele_freqs);
    N_alleles = size(summary.alleles,1);
    
    di_mu = effective_alleles(summary)/N_total;    
    
    datavec = repelem([1:N_alleles]', summary.allele_freqs);    
    
    di_trial = zeros(N_trials,1);
    for i = 1:N_trials
        sample_freqs = accumarray(datasample(datavec, N_total, 'Replace', true),1, [N_alleles, 1]);
        di_trial(i) = effective_alleles(sample_freqs, is_template)/N_total;
    end
    di_sig = sqrt(sum((di_trial-di_mu).^2)/N_trials);
end