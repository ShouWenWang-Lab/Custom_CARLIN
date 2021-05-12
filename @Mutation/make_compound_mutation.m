function m = make_compound_mutation(s, e, in_seq, ref_seq)

    in_seq = degap(in_seq);
    ref_seq = degap(ref_seq);
    min_len = min(length(in_seq), length(ref_seq));

    match_left = find(ref_seq(1:min_len)~=in_seq(1:min_len), 1, 'first')-1;
    if (isempty(match_left))
        match_left = min_len;
    end

    match_right = find(ref_seq(end:-1:end-min_len+1)~=in_seq(end:-1:end-min_len+1), 1, 'first')-1;
    if (isempty(match_right))
        match_right = min_len;
    end

    if (match_left+match_right < min_len)    
        max_len = max(length(in_seq), length(ref_seq));
        in_seq  = pad( in_seq(1+match_left:end-match_right), max_len-match_left-match_right, 'right', '-');
        ref_seq = pad(ref_seq(1+match_left:end-match_right), max_len-match_left-match_right, 'right', '-');
        m = Mutation('C', s+match_left, e-match_right, in_seq, ref_seq);
    else
        if (match_left+match_right > min_len)
            match_right = min_len-(match_left);
        end        
        if (match_left+match_right == length(in_seq))
            insert = ref_seq(1+match_left:end-match_right);
            m = Mutation('D', s+match_left, e-match_right, pad('', length(insert), '-'), insert);
        else
            if (match_left > 0)
                insert = in_seq(match_left:end-match_right);
                m = Mutation('I', s+match_left-1, e-match_right, insert, pad(ref_seq(match_left), length(insert), 'right', '-'));
            else
                insert = in_seq(1:end-match_right+1);
                m = Mutation('I', s, e-match_right+1, insert, pad(ref_seq(end-match_right+1), length(insert), 'left', '-'));
            end
        end
    end
end