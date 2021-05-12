classdef ExperimentSummary
    
    properties (SetAccess = immutable)
        alleles
        allele_freqs
    end
    
    methods (Access = public)
        
        function obj = ExperimentSummary(alleles, allele_freqs, preserve_order)
            
            % Always realign when we get to this stage. If we're merging datasets, 
            % makes sure results are updated, and makes sure that the same
            % sequence always provides the same allele regardless of
            % alignment.
            
            if (size(alleles,1) == 0)
                
                if (nargin > 1)
                    assert(size(allele_freqs,1)==0);
                end
                obj.alleles = {};
                obj.allele_freqs = [];
                
            else

                if (isa(alleles{1}, 'AlignedSEQ'))
                    alleles = cellfun(@(x) degap(x.get_seq), alleles, 'un', false);
                else
                    alleles = cellfun(@(x) degap(x), alleles, 'un', false);
                end

                if (nargin == 1)
                    [alleles, ~, obj.allele_freqs] = unique_by_freq(alleles);
                    obj.allele_freqs = accumarray(obj.allele_freqs,1);
                else
                    if (nargin == 2)
                        preserve_order = false;
                    end
                    assert(size(alleles,1) == size(allele_freqs,1));
                    if (issorted(allele_freqs,'descend') || preserve_order)
                        obj.allele_freqs = allele_freqs;
                    else
                        [obj.allele_freqs, idx] = sort(allele_freqs, 'descend');
                        alleles = alleles(idx);
                    end 
                end

                [~, obj.alleles] = cellfun(@(x) CARLIN_def.cas9_align(x), alleles, 'un', false);
            end
            
        end
        
    end
    
    methods (Static)
        
        function [obj, sample_map, allele_breakdown_by_sample] = FromMerge(samples)
            
            assert(isa(samples, 'ExperimentSummary'));
            if (size(samples,1)==1)                
                obj = ExperimentSummary(samples.alleles, samples.allele_freqs);
                sample_map = {[1:length(obj.alleles)]'};
                allele_breakdown_by_sample = samples.allele_freqs;
                return;
            end
            
            joint_frequency = vertcat(samples.allele_freqs);
            joint_alleles = {samples.alleles}';
            N_alleles = cellfun(@length, joint_alleles);            
            [joint_alleles, ~, sample_map] = ...
                unique_by_freq(cellfun(@(x) degap(x.get_seq), vertcat(joint_alleles{:}), 'un', false), joint_frequency);
            
            obj = ExperimentSummary(joint_alleles, accumarray(sample_map, joint_frequency));
            
            sample_map = mat2cell(sample_map, N_alleles);        
            N_samples = length(sample_map);
            allele_breakdown_by_sample = zeros(length(joint_alleles), N_samples);
            for i = 1:N_samples
                allele_breakdown_by_sample(sample_map{i},i) = samples(i).allele_freqs;
            end
        end
        
    end
end