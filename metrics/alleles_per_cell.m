function apc = alleles_per_cell(summary)
    assert(isa(summary, 'ExperimentSummary'));
    if (~isempty(summary.allele_freqs))
        apc = size(summary.alleles,1)/sum(summary.allele_freqs);
    else
        apc = 0;
    end
end
