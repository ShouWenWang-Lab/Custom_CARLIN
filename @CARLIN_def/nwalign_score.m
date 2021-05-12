function sc = nwalign_score(al, match_score, mismatch_penalty, gap_open, gap_extend)

    assert(isa(al, 'AlignedSEQ'));
    if (nargin < 5)
        gap_extend = 0.5;
    end
    if (nargin < 4)
        gap_open = 10;
    end
    if (nargin < 3)
        mismatch_penalty = -4;
    end
    if (nargin < 2)
        match_score = 5;
    end
    
    assert(mismatch_penalty <= 0 && all([match_score gap_open gap_extend] >= 0));
    
    seq = al.get_seq;
    ref = al.get_ref;
    
    seq_gap = isgap(seq);
    ref_gap = isgap(ref);
    
    assert(~any(seq_gap & ref_gap));
    
    both_nt = ~(seq_gap | ref_gap);
    
    num_match    = sum(seq(both_nt) == ref(both_nt));
    num_mismatch = sum(seq(both_nt) ~= ref(both_nt));
    
    del = (seq_gap & ~ref_gap);
    del = diff([0 del 0]);
    del = [find(del==1)', find(del==-1)'-1];
    ins = (~seq_gap & ref_gap);    
    ins = diff([0 ins 0]);
    ins = [find(ins==1)', find(ins==-1)'-1];
    indel = sortrows([del; ins]);
    assert(all(indel(2:end,1) > indel(1:end-1,2)));
    
    sc = num_match * match_score + num_mismatch * mismatch_penalty - sum(gap_open + gap_extend * diff(indel, [], 2));
    
end