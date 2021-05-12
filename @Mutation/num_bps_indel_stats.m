function [num_bp_del_mu, num_bp_ins_mu, num_bp_del_sd, num_bp_ins_sd] = num_bps_indel_stats(summary, edited_only)

    if (nargin < 2)
        edited_only = true;
    end

    assert(isa(summary, 'ExperimentSummary'));
    
    mut_events = cellfun(@Mutation.identify_Cas9_events, summary.alleles, 'un', false);
    [num_bp_del, num_bp_ins] = cellfun(@(x) arrayfun(@(i) x(i).num_bps_indel, [1:length(x)]'), mut_events, 'un', false);
    num_bp_del = cellfun(@(x) sum(x), num_bp_del);
    num_bp_ins = cellfun(@(x) sum(x), num_bp_ins);
    
    if (edited_only)
        
        is_template = ismember(cellfun(@(x) x.get_seq, summary.alleles, 'un', false), CARLIN_def.getInstance.seq.CARLIN);
        
        num_bp_ins_mu = sum(num_bp_ins(~is_template) .* summary.allele_freqs(~is_template))/sum(summary.allele_freqs(~is_template));
        num_bp_del_mu = sum(num_bp_del(~is_template) .* summary.allele_freqs(~is_template))/sum(summary.allele_freqs(~is_template));
        
        num_bp_ins_sd = sum(num_bp_ins(~is_template).^2 .* summary.allele_freqs(~is_template))/sum(summary.allele_freqs(~is_template));        
        
        num_bp_del_sd = sum(num_bp_del(~is_template).^2 .* summary.allele_freqs(~is_template))/sum(summary.allele_freqs(~is_template));        
        
    else
        num_bp_ins_mu = sum(num_bp_ins .* summary.allele_freqs)/sum(summary.allele_freqs);
        num_bp_del_mu = sum(num_bp_del .* summary.allele_freqs)/sum(summary.allele_freqs);
        
        num_bp_ins_sd = sum(num_bp_ins.^2 .* summary.allele_freqs)/sum(summary.allele_freqs);
        
        num_bp_del_sd = sum(num_bp_del.^2 .* summary.allele_freqs)/sum(summary.allele_freqs);
    end    
end