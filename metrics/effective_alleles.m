function ec = effective_alleles(varargin)
    if (nargin == 1)
        summary = varargin{1};
        assert(isa(summary, 'ExperimentSummary'));
        refseq = CARLIN_def.getInstance.seq.CARLIN;
        is_template = ismember(cellfun(@(x) degap(x.get_seq), summary.alleles, 'un', false), refseq);
        allele_counts = summary.allele_freqs;
    elseif (nargin == 2)
        allele_counts = varargin{1};
        is_template = varargin{2};
        assert(length(is_template) == length(allele_counts));
        assert(sum(is_template) <= 1);        
    end
    allele_counts = allele_counts(~is_template);
    allele_counts = allele_counts(allele_counts>0);
    p = allele_counts/sum(allele_counts);
    H = -sum(p.*log2(p));    
    ec = 2^H-isempty(p);    
end