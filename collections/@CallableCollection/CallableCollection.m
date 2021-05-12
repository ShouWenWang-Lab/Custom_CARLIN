classdef (Abstract) CallableCollection < handle
    
    methods (Abstract)
        call_result = call_alleles(obj, aligned, depth);
    end
    
    methods (Static)
        [alleles, which_seq, weight_contribution] = call_alleles_exact(aligned_seqs, aligned_seq_weights, top_only);
        [alleles, which_seq, weight_contribution] = call_alleles_coarse_grain(aligned_seqs, aligned_seq_weights, top_only);
    end
    
end