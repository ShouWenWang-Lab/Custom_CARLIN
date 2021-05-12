function [alleles, which_seqs, weight_contribution] = call_alleles_exact(aligned_seqs, aligned_seq_weights, dominant_only)

    if (nargin < 2)
        aligned_seq_weights = ones(size(aligned_seqs));
    else        
        assert(size(aligned_seqs,1) == size(aligned_seq_weights,1));
    end
    
    valid_mask = find(~cellfun(@isempty, aligned_seqs));
    aligned_seqs = aligned_seqs(valid_mask);
    aligned_seq_weights = aligned_seq_weights(valid_mask);    
    
    alleles = [];
    which_seqs = [];
    weight_contribution = [];
    
    if (isempty(aligned_seqs))
        return;
    end

    if (nargin < 3)
        dominant_only = true;
    end

    [seqs, ~, which_seq] = unique_by_freq(cellfun(@(x) x.get_seq(), aligned_seqs, 'un', false), aligned_seq_weights);
    seq_weight = accumarray(which_seq, aligned_seq_weights);
    assert(issorted(seq_weight, 'descend'));
    
    dominant_frac=0.5;
    if (dominant_only)
        if (seq_weight(1) / sum(seq_weight) >= dominant_frac)
            N = 1;
        else
            return;
        end
    else
        N = size(events,1);
    end
    
    alleles = cell(N,1);
    which_seqs = cell(N,1);
    weight_contribution = cell(N,1);
    
    for i = 1:N
        seq_mask = find(which_seq==i);
        [~, ref_ind] = max(aligned_seq_weights(seq_mask));
        ref_ind = seq_mask(ref_ind);
        alleles{i} = CARLIN_def.desemble_sequence(seqs{i}, aligned_seqs{ref_ind}.get_ref());
        which_seqs{i} = valid_mask(seq_mask);
        weight_contribution{i} = aligned_seq_weights(seq_mask);
    end
    
    if (N == 1)
        alleles = alleles{1};
        which_seqs = which_seqs{1};
        weight_contribution = weight_contribution{1};
    end
end