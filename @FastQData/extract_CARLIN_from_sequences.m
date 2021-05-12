function [SEQ_trimmed, read_SEQ_trimmed, masks, trim_loc] = extract_CARLIN_from_sequences(SEQ_raw, read_SEQ_raw, cfg)
            
    fprintf('Trimming CARLIN sequences from reads\n');
    
    ref = CARLIN_def.getInstance;
    
    % A one-off check to see if we match any CARLIN conserved sites exactly as a 
    % coarse way of quantifying off-target amplification. Not directly used in filtering.
    consite_loc = arrayfun(@(i) strfind(SEQ_raw, ref.getInstance.seq.consites{i}), 1:10, 'un', false);
    masks.CARLIN_match = any(~cellfun(@isempty, horzcat(consite_loc{:})),2);
    
    head_offset = 0;
    tail_offset = 0;
    if (startsWith(cfg.type, 'Bulk'))
        if (cfg.read_perspective.ShouldReverse == 'N')
            head_offset = cfg.UMI.length;
        else
            tail_offset = cfg.UMI.length;
        end
    end

    [masks.valid_5_primer, trim_loc.head_after_trim_5_primer] = ...
        FastQData.trim_at_scrutiny_level(cfg.trim.Primer5, SEQ_raw, ref.seq.Primer5, 'head', head_offset, ref.match_score.Primer5);

    [masks.valid_3_primer, trim_loc.tail_after_trim_3_primer] = ...
        FastQData.trim_at_scrutiny_level(cfg.trim.Primer3, SEQ_raw, ref.seq.Primer3, 'tail', tail_offset, ref.match_score.Primer3);

    [masks.valid_2_seq,    trim_loc.tail_after_trim_2_seq   ] = ...
        FastQData.trim_at_scrutiny_level(cfg.trim.SecondarySequence, SEQ_raw, ref.seq.SecondarySequence, ...
                                         'tail', ref.width.Primer3+tail_offset, ref.match_score.SecondarySequence);

    masks.valid_read_structure = ( masks.valid_5_primer & masks.valid_3_primer & masks.valid_2_seq );

    SEQ_trimmed = SEQ_raw;
    SEQ_trimmed(masks.valid_read_structure) = cellfun(@(x,s,e1,e2) x(s:min(e1,e2)), ...
                                                      SEQ_raw(masks.valid_read_structure), ...
                                                      num2cell(trim_loc.head_after_trim_5_primer(masks.valid_read_structure)), ...
                                                      num2cell(trim_loc.tail_after_trim_3_primer(masks.valid_read_structure)), ...
                                                      num2cell(trim_loc.tail_after_trim_2_seq(masks.valid_read_structure)), 'un', false);
    
    masks.trimmed_SEQ_long_enough = false(size(masks.valid_read_structure));
    masks.trimmed_SEQ_long_enough(masks.valid_read_structure) = ...
        cellfun(@length, SEQ_trimmed(masks.valid_read_structure)) >= ref.width.min_length;

    masks.SEQ_no_N = false(size(masks.valid_read_structure));
    masks.SEQ_no_N(masks.valid_read_structure) = ~cellfun(@(x) any(x=='N'), SEQ_trimmed(masks.valid_read_structure));

    masks.valid_SEQ_structure = masks.valid_read_structure & masks.trimmed_SEQ_long_enough & masks.SEQ_no_N;

    SEQ_trimmed(~masks.valid_SEQ_structure) = cellstr('');

    [SEQ_trimmed, ~, read_SEQ_trimmed] = unique_by_freq(SEQ_trimmed, accumarray(read_SEQ_raw,1));
    empty_loc = find(cellfun(@isempty, SEQ_trimmed));
    if (~isempty(empty_loc))
        SEQ_trimmed(empty_loc) = [];
        read_SEQ_trimmed(read_SEQ_trimmed == empty_loc) = 0;
        read_SEQ_trimmed(read_SEQ_trimmed > empty_loc) = read_SEQ_trimmed(read_SEQ_trimmed > empty_loc)-1;    
    end
    read_SEQ_trimmed = read_SEQ_trimmed(read_SEQ_raw);                                                        

    masks.CARLIN_match              = uint32(find(ismember(read_SEQ_raw, find(masks.CARLIN_match))));
    masks.valid_5_primer            = uint32(find(ismember(read_SEQ_raw, find(masks.valid_5_primer))));
    masks.valid_3_primer            = uint32(find(ismember(read_SEQ_raw, find(masks.valid_3_primer))));
    masks.valid_2_seq               = uint32(find(ismember(read_SEQ_raw, find(masks.valid_2_seq))));
    masks.valid_read_structure      = uint32(find(ismember(read_SEQ_raw, find(masks.valid_read_structure))));
    masks.trimmed_SEQ_long_enough   = uint32(find(ismember(read_SEQ_raw, find(masks.trimmed_SEQ_long_enough))));
    masks.SEQ_no_N                  = uint32(find(ismember(read_SEQ_raw, find(masks.SEQ_no_N))));
    masks.valid_SEQ_structure       = uint32(find(ismember(read_SEQ_raw, find(masks.valid_SEQ_structure))));

    N = length(read_SEQ_raw);

    fprintf('From %d reads, found (CARLIN,5,3,control,structure,long,AGCT,all) sequences (%d,%d,%d,%d,%d,%d,%d,%d) times\n', ...
        N, length(masks.CARLIN_match), length(masks.valid_5_primer), length(masks.valid_3_primer), length(masks.valid_2_seq), ...
        length(masks.valid_read_structure), length(masks.trimmed_SEQ_long_enough), length(masks.SEQ_no_N), length(masks.valid_SEQ_structure));  

end