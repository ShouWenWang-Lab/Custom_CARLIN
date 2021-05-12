function site = coarse_grain_motif(n)

    ref = CARLIN_def.getInstance;
    assert(n >= 1 && n <= ref.N.motifs);
    site = max(1,ceil((n-1)/3));
    
end