function [mut_list, bp_event, seq, ref_seq, ref_mask] = identify_sequence_events(aligned)

    % Naively identifies mutations just from relative alignment of sequence
    % and reference
    
    assert(isa(aligned, 'AlignedSEQ'));
    
    seq = aligned.get_seq();
    ref_seq = aligned.get_ref();
        
    ref = CARLIN_def.getInstance;    
    
    seq_gap = isgap(seq);
    ref_gap = isgap(ref_seq);
    assert(~any(seq_gap & ref_gap));
    
    ref_mask = find(~ref_gap);
    assert(length(ref_mask) == ref.width.CARLIN);
    which_bp = zeros(size(ref_seq));
    which_bp(ref_mask) = [1:ref.width.CARLIN];
    
    % 1. Default is no mutation    
    
    bp_event = repmat('N', [1, ref.width.CARLIN]);
    
    % 2. Tag deleted bps
    
    del = (seq_gap & ~ref_gap);
    del = diff([0 del 0]);
    del = [find(del==1)', find(del==-1)'-1];
    assert(all(arrayfun(@(b,e) all(ismembc([b:e]', ref_mask)), del(:,1), del(:,2))));
    
    del_mask = seq(ref_mask) == '-';
    bp_event(del_mask) = 'D';
    
    mut_list = {};
    if (~isempty(del))
        mut_list = arrayfun(@(bi,ei) Mutation('D', which_bp(bi), which_bp(ei), ...
                                               repmat('-', [1, ei-bi+1]), ref_seq(bi:ei)), del(:,1), del(:,2), 'un', false);
    end
    
    % 3. Tag inserted bps, special cases insertions that occur between the
    % cutsite and conserved sites, annotating them to be in the former.
    % I also include the bp from which the insertion happens in seq_old, 
    % aligned to the left or to the right according to the direction of the
    % insertion.
    
    frag_length = diff(ref_mask)-1;
    frag_ind = find(frag_length > 0);
    is = ismembc(frag_ind, ref.bounds.consites(:,2));
    bp_event(frag_ind(is)+1) = 'I';
    bp_event(frag_ind(~is))  = 'I';
    
    if (~isempty(frag_ind(is)))
        mut_list = [mut_list; arrayfun(@(i) Mutation('I', i+1, i+1, ...
                                                     seq(ref_mask(i)+1:ref_mask(i+1)), ...
                                                     [repmat('-', [1, frag_length(i)]) ref_seq(ref_mask(i+1))]), ...
                                                     frag_ind(is), 'un', false)'];
    end
    
    if (~isempty(frag_ind(~is)))
        mut_list = [mut_list; arrayfun(@(i) Mutation('I', i, i, ...                                       
                                                     seq(ref_mask(i):ref_mask(i+1)-1), ...
                                                     [ref_seq(ref_mask(i)) repmat('-', [1, frag_length(i)])]), ...
                                                     frag_ind(~is), 'un', false)'];
    end
        
    % Leading insertion
    if (ref_mask(1)~=1)
        bp_event(1) = 'I';
        mut_list = [mut_list; {Mutation('I', 1, 1, ...                                       
                               seq(1:ref_mask(1)), ...
                               [repmat('-', [1, ref_mask(1)-1]) ref_seq(ref_mask(1))])}];
    end
    
    % Trailing insertion
    if (ref_mask(end)~=length(seq))
        bp_event(end) = 'I';
        mut_list = [mut_list; {Mutation('I', ref.width.CARLIN, ref.width.CARLIN, ...                                       
                                        seq(ref_mask(end):end), ...
                                        [ref_seq(ref_mask(end)) repmat('-', [1, length(seq)-ref_mask(end)])])}];
    end
    
    % 4. Tag substituted bps
        
    bp_event(seq(ref_mask)~=ref_seq(ref_mask) & bp_event=='N') = 'M';    
    mut_list = [mut_list; arrayfun(@(i) Mutation('M', i, i, seq(ref_mask(i)), ref.seq.CARLIN(i)), find(bp_event=='M'), 'un', false)'];
    
    if (~isempty(mut_list))
        mut_list = vertcat(mut_list{:});
        [~, idx] = sortrows([vertcat(mut_list.loc_start), vertcat(mut_list.loc_end)], [1, 2]);
        mut_list = mut_list(idx);
        assert(all(vertcat(mut_list.loc_start) <= vertcat(mut_list.loc_end)));
    end
end