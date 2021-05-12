classdef UMIData < CallableCollection
    
    properties(SetAccess = immutable)
        UMI
        SEQ_ind
        SEQ_weight
    end
    
    methods (Access = public)
        
        function obj = UMIData(UMI, SEQ_ind, SEQ_weight)
            assert(size(SEQ_ind,1) == size(SEQ_weight,1));  
            assert(all(SEQ_weight > 0));
            assert(length(SEQ_ind) == length(unique(SEQ_ind)));
            obj.UMI = UMI;
            obj.SEQ_ind = SEQ_ind;
            obj.SEQ_weight = SEQ_weight;
        end
        
        function call_result = call_alleles(obj, depth, aligned)
            
            assert(isscalar(depth) && depth > 0);
            call_result.allele = [];
            call_result.constituents = [];
            
            if (sum(obj.SEQ_weight) >= depth)
                
                assert(isa(aligned, 'AlignedSEQDepot'));                    
                aligned_SEQ = arrayfun(@(x) aligned.get_alignment_for_SEQ_ind(x), obj.SEQ_ind, 'un', false);                                
                [call_result.allele, call_result.constituents, weight_contribution] ...
                    = CallableCollection.call_alleles_coarse_grain(aligned_SEQ, obj.SEQ_weight, true);
                
                if (~isempty(call_result.allele))
                    
                    assert(sum(obj.SEQ_weight(call_result.constituents)) == sum(weight_contribution));
                    event = call_result.allele.get_event_structure;
                    
                    if (startsWith(event, 'E') || endsWith(event, 'E'))
                        call_result.allele = [];
                        call_result.constituents = [];
                    end                    
                end
            end
        end        
    end
    
    methods (Static)
        
        function umi_data = FromFQ(UMI, FQ, CB_filter)
            assert(isa(FQ, 'FastQData'));
            if (isa(FQ, 'BulkFastQData'))
                assert(nargin == 2);                
                [SEQ_ind, SEQ_weight] = FQ.get_SEQ_ind_by_UMI(UMI);                
            elseif (isa(FQ, 'SCFastQData'))
                assert(nargin == 3);                
                [SEQ_ind, SEQ_weight] = FQ.get_SEQ_ind_by_UMI(UMI, CB_filter);
            end
            umi_data = UMIData(UMI, SEQ_ind, SEQ_weight);
        end
        
        function umi_data = FromMerge(UMI, UMIs_to_merge)
            merged_SEQ_weight = vertcat(UMIs_to_merge.SEQ_weight);
            [merged_SEQ_ind, ~, ind] = unique(vertcat(UMIs_to_merge.SEQ_ind));
            merged_SEQ_weight = accumarray(ind, merged_SEQ_weight);
            umi_data = UMIData(UMI, merged_SEQ_ind, merged_SEQ_weight);
        end 
        
    end
end
        