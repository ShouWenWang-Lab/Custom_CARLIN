function chop_sites = partition_sequence(seq_aligned, ref_aligned)

    ref = CARLIN_def.getInstance;

    assert(length(seq_aligned) == length(ref_aligned));
    assert(sum(ref_aligned~='-') == ref.width.CARLIN);
    
    ref_sites = ref.bounds.ordered;

    bp_pos = find(ref_aligned ~= '-');
    chop_sites = bp_pos(ref_sites);

    chop_sites(1,1) = 1;
    chop_sites(end,2) = length(ref_aligned);

    for i=1:size(chop_sites,1)-1
        if (chop_sites(i,2) < chop_sites(i+1,1)-1)
            % Gap between prefix/PAM and consite. Can assign to either
            % arbitrarily. Pick consite as convention.
            if (ismember(i, ref.motifs.prefix) || ismember(i, ref.motifs.pams))
                chop_sites(i+1,1) = chop_sites(i,2)+1;
            % Gap between consite and cutsite. Merge gap into cutsite
            elseif (ismember(i, ref.motifs.consites))
                chop_sites(i+1,1) = chop_sites(i,2)+1;
            % Gap between cutsite and PAM/postfix.
            % Merge gap into cutsite
            elseif (ismember(i, ref.motifs.cutsites))
                chop_sites(i,2) = chop_sites(i+1,1)-1;
            end
        end
    end
    
end