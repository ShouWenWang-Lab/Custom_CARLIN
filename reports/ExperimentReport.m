classdef (Abstract=true) ExperimentReport < ExperimentSummary
    
    properties (SetAccess = immutable)
        reads
        allele_colony
        N
    end
    
    methods (Access = public)
        
        function obj = ExperimentReport(alleles, allele_colony, N, reads)
            assert(size(alleles,1) == size(allele_colony,1));
            obj = obj@ExperimentSummary(alleles, cellfun(@length, allele_colony));
            obj.allele_colony = allele_colony;
            obj.N = N;
            obj.reads = reads;
        end
        
        function conv = ExperimentSummary(obj)
           conv = ExperimentSummary(obj.alleles, obj.allele_freqs, true); 
        end
        
    end
    
end