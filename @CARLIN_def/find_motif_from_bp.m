function motif = find_motif_from_bp(n)
    
    ref = CARLIN_def.getInstance;
    assert(n >= 1 && n <= ref.width.CARLIN);
    motif = find(ref.bounds.ordered(:,1) <= n & ref.bounds.ordered(:,2) >= n);
    
end