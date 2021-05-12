classdef (Sealed=true) BulkFastQData < FastQData
    
    methods (Static)
        
        [SEQ, read_SEQ, QC, Nreads] = parse_bulk_fastq(fastq_file, cfg);
        [UMI, read_UMI, SEQ, read_SEQ, QC, trim_loc] = parse_bulk_lines(SEQ, read_SEQ, QC, cfg, trim_loc);
        masks = filter_bulk_UMIs(UMI, read_UMI, QC, L,QC_threshold) % SW: I add 'QC_threshold' here            
         
    end 
        
    methods (Access = public)
    
        function obj = BulkFastQData(fastq_file, cfg)
            
            assert(nargin == 2, 'Expected two inputs to BulkFastQData constructor');
            
            % Logical flow is a little convoluted because UMIs and CARLIN are part of the 
            % same read and we still want to get good diagnostics out
            
            % 1. Just get reads (SEQ=UMI+CARLIN), and QC
            
            [SEQ, read_SEQ, QC, Nreads] = BulkFastQData.parse_bulk_fastq(fastq_file, cfg);
            
            % 2. Trim out CARLIN based on Primer5 and Primer3 settings. 
            
            [SEQ_trimmed, read_SEQ_trimmed, seq_masks, trim_loc] = FastQData.extract_CARLIN_from_sequences(SEQ, read_SEQ, cfg);
            
            % 3. What's left over after trimming (specified by trim_loc) is
            % used to get UMIs. Reconstitute 'raw' SEQ by keeping what's
            % left of the original read after excising UMIs. trim_loc
            % indexes this new SEQ_raw after trimming for UMIs
            
            [UMI, read_UMI, SEQ, read_SEQ, QC, trim_loc] = BulkFastQData.parse_bulk_lines(SEQ, read_SEQ, QC, cfg, trim_loc);
            
            header_masks = BulkFastQData.filter_bulk_UMIs(UMI, read_UMI, QC, cfg.UMI.length,cfg.UMI.QC);
                                               
            masks = merge_structs(header_masks, seq_masks);
            masks.valid_lines = intersect(masks.valid_provenance_structure, masks.valid_SEQ_structure);
            fprintf('Merging filters\n');
            fprintf('From %d reads, found valid (provenance,sequence,both) reads (%d,%d,%d) times\n', ...
                Nreads, length(masks.valid_provenance_structure), length(masks.valid_SEQ_structure), length(masks.valid_lines));
                        
            obj = obj@FastQData(fastq_file, Nreads, SEQ, read_SEQ, SEQ_trimmed, read_SEQ_trimmed, UMI, read_UMI, QC, masks, trim_loc);
            
        end
    end
end
