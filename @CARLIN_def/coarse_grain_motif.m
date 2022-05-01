function site = coarse_grain_motif(n)

    ref = CARLIN_def.getInstance;
    assert(n >= 1 && n <= ref.N.motifs);
    site = max(1,ceil((n-1)/3));
    
    % 20220430: fix the upper bound manually. This deals with the postfix motif
    % this change is introduced to account for RC and TC, which has one more pam motif
    % that that for the def_CARLIN_cCARLIN
    if site>10
        site=10;
    end
end