function spc = singletons_per_cell(summary)
    assert(isa(summary, 'ExperimentSummary'));
    if (~isempty(summary.allele_freqs))
        spc = sum(summary.allele_freqs==1)/sum(summary.allele_freqs);
    else
        spc = 0;
    end
end