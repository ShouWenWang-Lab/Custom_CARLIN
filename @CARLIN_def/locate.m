function loc = locate(ref, bp, motif_bounds)

    loc.abs = find(motif_bounds(:,1) <= bp & motif_bounds(:,2) >= bp);    
    assert(~isempty(loc.abs) && loc.abs >= 1 && loc.abs <= ref.N.motifs);
    loc.pos = bp-motif_bounds(loc.abs,1)+1;
    loc.rel = round(loc.abs/3);
    
    if (loc.abs==1)
        loc.type = 'prefix';
        loc.rel = 1;
    elseif (loc.abs==ref.N.motifs)
        loc.type = 'postfix';
        loc.rel = 1;
    elseif (mod(loc.abs,3)==0)
        loc.type = 'cutsites';
    elseif (mod(loc.abs,3)==2)
        loc.type = 'consites';
    elseif (mod(loc.abs,3)==1)
        loc.type = 'pams';
    end 
end