classdef BulkUMICollection < TaggedCollection & CallableCollection
    
    properties(SetAccess = immutable)
        UMIs
    end
    
    methods (Access = public)
        
        function obj = BulkUMICollection(UMI_data)
            obj.UMIs = UMI_data;
        end
        
        function call_result = call_alleles(obj, depth, aligned)
            
            assert(isscalar(depth) && depth > 0);
            
            if (depth > 1)
                fprintf('Calling alleles for UMI collection\n');
            end
            
            assert(isa(aligned, 'AlignedSEQDepot'));
            N = size(obj.UMIs,1);
            call_result = cell(N,1);
            temp = obj.UMIs;
            
            % This can be parfor on a bigmem machine
            for i = 1:N            
                call_result{i} = temp(i).call_alleles(depth, aligned);
            end
            
            uncalled = cellfun(@isempty, call_result);
            call_result = vertcat(call_result{:});
            temp = {obj.UMIs(~uncalled).UMI};
            [call_result.UMI] = temp{:};
            call_result(cellfun(@isempty, {call_result.allele}')) = [];
            
            if (depth > 1)
	            fprintf('...%d UMIs callable\n', length(call_result));
            end
            
        end
        
        function [out, UMI_map] = denoise(obj, aligned)
            
            fprintf('Denoising UMI collection\n');
            k = {obj.UMIs.UMI}';
            v = k;
            
            call_result = obj.call_alleles(1, aligned);
            [is, where] = ismember(k, {call_result.UMI}');
            valid_UMIs = find(is);
            
            weight_for_UMIs = cellfun(@(x) sum(x), {obj.UMIs(valid_UMIs).SEQ_weight}');
            
            [alleles, ~, allele_for_UMIs] = unique(cellfun(@(x) x.get_seq(), {call_result(where(is)).allele}', 'un', false));
            
            N = size(alleles,1);
            allele_mask = cell(N,1);
            UMIs_for_allele = cell(N,1);
            UMI_map = cell(N,1);            

            parfor i = 1:N
                allele_mask{i} = allele_for_UMIs==i;
                UMIs_for_allele{i} = k(valid_UMIs(allele_mask{i}));
                UMI_map{i} = zeros(size(UMIs_for_allele{i}));
                weight_for_allele = weight_for_UMIs(allele_mask{i});
                UMI_map{i} = TaggedCollection.directional_adjacency_top_down_denoiser(UMIs_for_allele{i}, weight_for_allele);
                assert(all(weight_for_allele <= weight_for_allele(UMI_map{i})));
            end
            
            for i = 1:N
                v(valid_UMIs(allele_mask{i})) = UMIs_for_allele{i}(UMI_map{i});
            end
            
            UMI_map = containers.Map(k, v);
            [v, ~, ind] = unique(v);
            N = size(v,1);
            UMI_data= cell(N,1);
            merged_data = arrayfun(@(i) vertcat(obj.UMIs(ind==i)), 1:N, 'un', false);

            parfor i = 1:N
                UMI_data{i} = UMIData.FromMerge(v{i}, merged_data{i});
            end
            out = BulkUMICollection(vertcat(UMI_data{:}));
            fprintf('...from %d UMIs to %d UMIs\n', length(obj.UMIs), N);
        end
        
        function UMI_freq = get_UMI_freq(obj)
            UMI_freq = sort(cellfun(@(x) sum(x), {obj.UMIs.SEQ_weight}'), 'descend');
        end
        
        function thresholds = compute_thresholds(obj, params, FQ)
            fprintf('Computing thresholds for UMI collection\n');            
            p = cellfun(@(x) 10.^((double(x)-33)/-10), FQ.QC(FQ.masks.valid_lines), 'un', false);
            p = sort(horzcat(p{:}), 'descend');
            p = p(ceil(length(p)/20));
            UMI_freq = obj.get_UMI_freq();
            thresholds = TaggedCollection.threshold_function(UMI_freq, params.Results.max_molecules, ...
                                                             length(FQ.masks.valid_lines), p, ...
                                                             params.Results.read_cutoff_UMI_denoised, ...
                                                             params.Results.read_override_UMI_denoised, 'UMI');
            fprintf('...%d common UMIs\n', sum(UMI_freq >= thresholds.chosen));
        end
    
    end
    
    methods (Static)
        
        function out = FromFQ(FQ)            
            assert(isa(FQ, 'BulkFastQData'));
            UMIs = FQ.get_UMIs();            
            N = size(UMIs,1);
            fprintf('Building UMI collection from FastQ with %d UMIs\n', N);
            [SEQ_ind, SEQ_weight] = FQ.get_SEQ_ind_by_UMI(UMIs);
            assert(size(SEQ_ind,1) == N && size(SEQ_weight,1) == N);            
            UMI_data = cellfun(@(u,i,w) UMIData(u,i,w), UMIs, SEQ_ind, SEQ_weight, 'un', false);            
            out = BulkUMICollection(vertcat(UMI_data{:}));
        end 
        
    end
end