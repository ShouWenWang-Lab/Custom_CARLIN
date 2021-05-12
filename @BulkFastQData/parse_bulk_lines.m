function [UMI, read_UMI, SEQ, read_SEQ, QC, trim_loc] = parse_bulk_lines(SEQ, read_SEQ, QC, cfg, trim_loc)
    
    assert(size(QC,1) == size(read_SEQ,1), 'Length of QC and SEQs must match when trimming for UMIs');    
    assert(all(structfun(@(x) size(x,1) == size(SEQ,1), trim_loc)), 'Length of trim_loc and SEQs must match when trimming for UMIs');
    
    ref = CARLIN_def.getInstance;
    
    if (cfg.read_perspective.ShouldReverse=='N')
        end_loc = trim_loc.head_after_trim_5_primer-length(ref.seq.Primer5);
        start_loc = end_loc-cfg.UMI.length;
        start_loc = max(1, start_loc);
        end_loc = max(0, end_loc-1);
        trim_loc = structfun(@(x) x-end_loc, trim_loc, 'un', false);
        [UMI, SEQ] = cellfun(@(x,s,e) deal(x(s:e), x(e+1:end)), SEQ, num2cell(start_loc), num2cell(end_loc), 'un', false);
    else
        start_loc = trim_loc.tail_after_trim_3_primer+length(ref.seq.Primer3);
        end_loc = start_loc+cfg.UMI.length;
        start_loc = min(cellfun(@length, SEQ)+1,double(start_loc)+1);
        end_loc = min(cellfun(@length, SEQ), double(end_loc));
        [SEQ, UMI] = cellfun(@(x,s,e) deal(x(1:s-1), x(s:e)), SEQ, num2cell(start_loc), num2cell(end_loc), 'un', false);
        
    end
    
    QC = cellfun(@(x,s,e) x(s:e), QC, num2cell(start_loc(read_SEQ)), num2cell(end_loc(read_SEQ)), 'un', false);
    
    [UMI, QC] = FastQData.orient_reads(cfg, UMI, QC);
    
    [UMI, ~, ind_UMI] = unique_by_freq(UMI, accumarray(read_SEQ,1));
    
    [SEQ, backtrack, ind_SEQ] = unique_by_freq(SEQ, accumarray(read_SEQ,1));
    trim_loc = structfun(@(x) x(backtrack), trim_loc, 'un', false);    
    
    read_UMI = ind_UMI(read_SEQ);
    read_SEQ = ind_SEQ(read_SEQ);
       
end