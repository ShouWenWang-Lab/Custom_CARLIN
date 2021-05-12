classdef SCExperimentReport < ExperimentReport
    
    methods (Access = public)
        
        function obj = SCExperimentReport(CB_collection, CB_collection_denoised, CB_denoise_map, CB_called_allele, FQ, threshold, ref_CBs)
            
            fprintf('Generating report for SC experiment\n');
            
            [alleles, ~, which_allele] = unique_by_freq(cellfun(@(x) degap(x.get_seq), {CB_called_allele.allele}', 'un', false));            
            allele_colony = arrayfun(@(i) {CB_called_allele(which_allele==i).CB}', [1:size(alleles,1)]', 'un', false);
            
            template_index = strcmp(alleles, CARLIN_def.getInstance.seq.CARLIN);
            temp = {CB_collection.CBs.UMIs}';
            cb_weights = zeros(size(temp));
            cb_weights(~cellfun(@isempty, temp)) = cellfun(@(x) sum(vertcat(x.SEQ_weight)), temp(~cellfun(@isempty, temp)));
            
            N.reference_tags = length(ref_CBs);
            N.uncleaned_tags = length(CB_collection_denoised.CBs);
            N.matched_tags   = length(keys(CB_denoise_map.CB));            
            N.cleaned_tags   = length(unique(values(CB_denoise_map.CB)));
            N.common_tags    = sum(cb_weights >= threshold.CB.chosen);
            N.called_tags    = size(CB_called_allele,1);
            N.eventful_tags  = sum(cellfun(@length, allele_colony(~template_index)));
                        
            reads.in_fastq            = FQ.Nreads;
            for fn = fieldnames(FQ.masks)'
                reads.(fn{1}) = length(FQ.masks.(fn{1}));
            end
            [is, where] = ismember(keys(CB_denoise_map.CB), FQ.CB);
            assert(all(is));
            reads.matched_tags         = sum(ismember(FQ.read_CB(FQ.masks.valid_lines), where));
            reads.common_tags          = sum(cb_weights(cb_weights >= threshold.CB.chosen));
            
            [is, where] = ismember({CB_called_allele.CB}', {CB_collection.CBs.CB}');
            assert(all(is));
            reads.called_tags_total  = cb_weights(where);
            reads.called_tags_allele = zeros(length(where),1);
            
            for i = 1:length(where)
                reads.called_tags_allele(i) = sum(arrayfun(@(j) sum(CB_collection.CBs(where(i)).UMIs(j).SEQ_weight ...
                                                  (CB_called_allele(i).umi_call_result(j).constituents)), CB_called_allele(i).constituents));
            end
            
            if (any(template_index))
                is = ismember({CB_called_allele.CB}', allele_colony{template_index});
                reads.eventful_tags_total  = sum(reads.called_tags_total(~is));
                reads.eventful_tags_allele = sum(reads.called_tags_allele(~is));
            else
                reads.eventful_tags_total = sum(reads.called_tags_total);
                reads.eventful_tags_allele = sum(reads.called_tags_allele);
            end
            
            reads.called_tags_total = sum(reads.called_tags_total);
            reads.called_tags_allele = sum(reads.called_tags_allele);
            
            obj = obj@ExperimentReport(alleles, allele_colony, N, reads);
            
        end
        
    end
end