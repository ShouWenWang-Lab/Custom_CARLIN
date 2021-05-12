classdef CBCollection < TaggedCollection & CallableCollection
    
    properties(SetAccess = immutable)
        CBs
    end
    
    methods (Access = public)
        
        function obj = CBCollection(CB_data)
            obj.CBs = CB_data;
        end
        
        function call_result = call_alleles(obj, depth, aligned)
            
            fprintf('Calling alleles for CB collection\n');
            assert(isequal(size(depth), [1, 2]) && all(depth > 0) && (depth(1) >= depth(2)));                        
            assert(isa(aligned, 'AlignedSEQDepot'));
            
            N = size(obj.CBs,1);
            call_result = cell(N,1);
            temp = obj.CBs;
            
            % This can be parfor on a bigmem machine
            for i = 1:N
                call_result{i} = temp(i).call_alleles(depth, aligned);
            end
            
            % Call result has a per-CB allele (which may or may not be empty) and a 
            % per-UMI call result (which may or may not be
            % empty, even when per-CB allele is empty).
            
            % Get cell array elements where all fields are empty.
            uncalled = cellfun(@isempty, call_result);
            
            % This concatenation collapses these all-empty elements
            call_result = vertcat(call_result{:});
            
            if (isempty(call_result))
                temp = {'CB', 'allele'};
                temp{2,1} = {};
                call_result = struct(temp{:});
            else
                % Get CBs of uncollapsed elements and tack them on as a field
                temp = {obj.CBs(~uncalled).CB};
                [call_result.CB] = temp{:};

                % Remove CBs where there is no callable per-CB allele
                call_result(cellfun(@isempty, {call_result.allele}')) = [];

                % Corner case, if nothing is callable, make me a 0x1 array
                % instead of a 1x0 array
                if (isempty(call_result))
                    call_result = call_result';
                end
            end
            fprintf('...%d CBs callable\n', size(call_result,1));
        end
        
        function [out, collapse_map] = denoise(obj, ref_list)
            
            fprintf('Denoising CB collection\n');

            [is, which_ref] = ismember({obj.CBs.CB}', ref_list);
            CB_weight = cellfun(@(x) sum(vertcat(x.SEQ_weight)), {obj.CBs.UMIs}');
            CB_MiSEQ_internal = TaggedCollection.directional_adjacency_top_down_denoiser({obj.CBs.CB}', CB_weight, is);
            
            self_same = CB_MiSEQ_internal==[1:length(CB_MiSEQ_internal)]';

            matched = ismember(CB_MiSEQ_internal,find(self_same&is));
            CB_cleaned = cell(size(CB_weight));
            CB_cleaned(matched) = ref_list(which_ref(CB_MiSEQ_internal(matched)));

            matched = find(matched);
            k = {obj.CBs(matched).CB};
            v = CB_cleaned(matched);
            collapse_map.CB = containers.Map(k, v);
            [v, ~, ind] = unique(v);
            N = size(v,1);
            
            fprintf('...(%d/%d) MiSEQ CBs matched to (%d/%d) reference CBs\n', ...
                length(matched), length(obj.CBs), N, length(ref_list));
            
            fprintf('Denoising UMIs\n')
            CB_data = cell(N,1);
            UMI_map = cell(N,1); 
            merged_data = arrayfun(@(i) vertcat(obj.CBs(matched(ind==i))), 1:N, 'un', false);            
            
            parfor i = 1:N
                temp = CBData.FromMerge(v{i}, merged_data{i});
                [CB_data{i}, UMI_map{i}] = temp.denoise();
            end
            collapse_map.UMI = containers.Map(v, UMI_map);
            out = CBCollection(vertcat(CB_data{:}));
            
        end
        
        function CB_freq = get_CB_freq(obj)
            CB_freq = sort(cellfun(@(x) sum(vertcat(x.SEQ_weight)), {obj.CBs.UMIs}'), 'descend');
        end
        
        function UMI_freq = get_UMI_freq(obj)
            UMI_freq = cellfun(@(x) cellfun(@(y) sum(y), {x.SEQ_weight}'), {obj.CBs.UMIs}', 'un', false);
            UMI_freq = sort(vertcat(UMI_freq{:}), 'descend');
        end
        
        function thresholds = compute_thresholds(obj, params, FQ, N_ref_CBs)
            
            fprintf('Computing thresholds for CB collection\n');
            p = cellfun(@(x) 10.^((double(x)-33)/-10), FQ.QC(FQ.masks.valid_lines), 'un', false);
            p = sort(horzcat(p{:}), 'descend');
            p = p(ceil(length(p)/20));
            CB_freq = obj.get_CB_freq();
            thresholds.CB = TaggedCollection.threshold_function(CB_freq, min(params.Results.max_cells, N_ref_CBs), ...
                                                                length(FQ.masks.valid_lines), p, ...
                                                                params.Results.read_cutoff_CB_denoised, ...
                                                                params.Results.read_override_CB_denoised, 'CB');
            fprintf('...%d common CBs\n', sum(CB_freq >= thresholds.CB.chosen));
            
            UMI_freq = obj.get_UMI_freq();
            thresholds.UMI = TaggedCollection.threshold_function(UMI_freq, inf, ...
                                                                 length(FQ.masks.valid_lines), p, ...
                                                                 params.Results.read_cutoff_UMI_denoised, ...
                                                                 params.Results.read_override_UMI_denoised, 'UMI');
            fprintf('...%d common UMIs\n', sum(UMI_freq >= thresholds.UMI.chosen));
                                           
        end
        
    end
    
    methods (Static)
        
        function out = FromFQ(FQ)
            
            assert(isa(FQ, 'SCFastQData'))
            CBs = FQ.get_CBs();
            N = size(CBs,1);
            fprintf('Building CB collection from FastQ with %d CBs\n', N);
            
            [UMIs, filter] = FQ.get_UMIs_by_CB(CBs);
            CB_data = cell(N,1);
            
            % This can be parfor on a bigmem machine
            for i = 1:N
                [SEQ_ind, SEQ_weight] = FQ.get_SEQ_ind_by_UMI(UMIs{i}, filter{i});
                UMI_data = cellfun(@(u,i,w) UMIData(u,i,w), UMIs{i}, SEQ_ind, SEQ_weight, 'un', false);
                CB_data{i} = CBData(CBs{i}, vertcat(UMI_data{:}));
            end
            out = CBCollection(vertcat(CB_data{:}));
        end
        
    end
end