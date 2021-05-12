classdef (Sealed=true) SCFastQData < FastQData
    
    properties (SetAccess = immutable, GetAccess = public)
        CB
        read_CB
    end
    
    methods (Static)
        
        masks = filter_sc_CBs_and_UMIs(cfg, CB, read_CB, UMI, read_UMI, QC);
        [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = parse_indrops_fastq(fastq_file);
        [CB, UMI, QC] = parse_indrops_provenance(H)
        [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = parse_10x_fastq(fastq_file, cfg);
        [CB, UMI, QC] = parse_10x_provenance(CB, QC, cfg);
        
    end
        
    methods (Access = public)
        
        function obj = SCFastQData(fastq_file, cfg)
            
            assert(nargin >= 2, 'Expected two inputs to SCFastQDat constructor');
            
            assert(strcmp(cfg.type, 'SC'), 'Invalid CFG passed function');
            if(strcmp(cfg.SC.Platform, 'InDrops'))
                [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = SCFastQData.parse_indrops_fastq(fastq_file);
            elseif(strcmp(cfg.SC.Platform, '10x'))
                [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = SCFastQData.parse_10x_fastq(fastq_file, cfg);            
            else
                error('Unsupported SC Platform: %s', cfg.Platform);
            end
            
            header_masks = SCFastQData.filter_sc_CBs_and_UMIs(cfg, CB, read_CB, UMI, read_UMI, QC);
            [SEQ_trimmed, read_SEQ_trimmed, seq_masks, trim_loc] = FastQData.extract_CARLIN_from_sequences(SEQ, read_SEQ, cfg);
            
            masks = merge_structs(header_masks, seq_masks);
            masks.valid_lines = intersect(header_masks.valid_provenance_structure, seq_masks.valid_SEQ_structure);
            fprintf('Merging filters\n');
            fprintf('From %d reads, found valid (provenance,sequence,both) reads (%d,%d,%d) times\n', ...
                Nreads, length(masks.valid_provenance_structure), length(masks.valid_SEQ_structure), length(masks.valid_lines));
            
            obj = obj@FastQData(fastq_file, Nreads, SEQ, read_SEQ, SEQ_trimmed, read_SEQ_trimmed, UMI, read_UMI, QC, masks, trim_loc);
            obj.CB = CB;
            obj.read_CB = uint32(read_CB);
        end
        
        function [UMIs, filter] = get_UMIs_by_CB(obj, CB)
            
            N = size(CB,1);
            if (N > 1)                
                assert(iscell(CB));
            else
                CB = cell(CB);
            end
            
            filter = cell(N,1);
            UMIs = cell(N,1);
            [is, which_CB] = ismember(CB, obj.CB);
            
            if (any(is))
                temp = obj.read_CB(obj.masks.valid_lines);
                filter(is) = arrayfun(@(i) obj.masks.valid_lines(temp==i), which_CB(is), 'un', false);
                UMIs(is) = cellfun(@(x) obj.UMI(unique_by_freq(obj.read_UMI(x))), filter(is), 'un', false);
            end
            
            filter(~is) = {[]};
            UMIs(~is) = {[]};
        end
        
        function CBs = get_CBs(obj)
            CBs = obj.CB(unique_by_freq(obj.read_CB(obj.masks.valid_lines)));
        end
        
    end
end
