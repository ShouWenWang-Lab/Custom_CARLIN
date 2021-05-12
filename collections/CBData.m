classdef CBData < TaggedCollection & CallableCollection
 
    properties(SetAccess = immutable)
        CB
        UMIs
    end        

    methods (Access = public)
        
         function obj = CBData(CB, UMI_data)             
             obj.CB = CB;
             obj.UMIs = UMI_data;
         end
        
        function call_result = call_alleles(obj, depth, aligned)
            
            assert(isequal(size(depth), [1, 2]) && all(depth > 0) && (depth(1) >= depth(2)));
            call_result = [];
            
            if (~isempty(obj.UMIs))
                umi_weights = vertcat(obj.UMIs.SEQ_weight);
                if (sum(umi_weights) >= depth(1))
                    assert(isa(aligned, 'AlignedSEQDepot'));
                    call_result.umi_call_result = cell(size(obj.UMIs,1),1);
                    for i = 1:size(obj.UMIs,1)
                        call_result.umi_call_result{i} = obj.UMIs(i).call_alleles(depth(2), aligned);
                    end
                    call_result.umi_call_result = vertcat(call_result.umi_call_result{:});
                    [call_result.umi_call_result.UMI] = obj.UMIs.UMI;
                    [call_result.allele, call_result.constituents] = ...
                        CallableCollection.call_alleles_coarse_grain({call_result.umi_call_result(:).allele}');
                end
            end
        end
        
        function [out, UMI_map] = denoise(obj)
            
            if (isempty(obj.UMIs))
                out = CBData(obj.CB, []);
                UMI_map = containers.Map();
                return;
            end
            
            k = {obj.UMIs.UMI}';                      
            weight_for_UMIs = arrayfun(@(i) sum(obj.UMIs(i).SEQ_weight), [1:size(k,1)]');
            UMI_map = TaggedCollection.directional_adjacency_top_down_denoiser(k, weight_for_UMIs);
            assert(all(weight_for_UMIs <= weight_for_UMIs(UMI_map)));
            v = k(UMI_map);

            UMI_map = containers.Map(k, v);
            [v, ~, ind] = unique(v);
            UMI_data = cell(size(v,1),1);
            for i = 1:size(v,1)
                UMI_data{i} = UMIData.FromMerge(v{i}, vertcat(obj.UMIs(ind==i)));
            end
            out = CBData(obj.CB, vertcat(UMI_data{:}));
        end
        
    end
    
    methods (Static)
        
        function out = FromFQ(CB, FQ)
            assert(isa(FQ, 'SCFastQData'));
            [UMIs, lines] = FQ.get_UMIs_by_CB(CB);            
            N = size(UMIs,1);
            UMI_data = cell(N,1);
            for i = 1:size(UMIs,1)
                UMI_data{i} = UMIData.FromFQ(UMIs{i}, FQ, lines{i});
            end
            out = CBData(CB, vertcat(UMI_data{:}));
        end
        
        function CB_data = FromMerge(CB, data_to_merge)            
            data_to_merge = vertcat(data_to_merge.UMIs);
            if (isempty(data_to_merge))
                CB_data = CBData(CB, []);
            else
                [merged_UMIs, ~, ind] = unique(cellstr(vertcat(data_to_merge.UMI)));
                N = size(merged_UMIs,1);
                UMI_data = cell(N,1);
                for i = 1:N
                    UMI_data{i} = UMIData.FromMerge(merged_UMIs{i}, data_to_merge(ind==i));
                end
                CB_data = CBData(CB, vertcat(UMI_data{:}));
            end
        end
        
    end
end