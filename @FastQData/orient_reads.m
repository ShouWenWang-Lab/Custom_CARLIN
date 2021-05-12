function [SEQ, QC] = orient_reads(cfg, SEQ, QC)

    if (cfg.read_perspective.ShouldReverse=='Y')
        SEQ = cellfun(@fliplr, SEQ, 'UniformOutput', false);
        QC  = cellfun(@fliplr, QC, 'UniformOutput', false);
    end
    
    if (cfg.read_perspective.ShouldComplement=='Y')
        SEQ = cellfun(@seqcomplement, SEQ, 'UniformOutput', false);
    end
    
end