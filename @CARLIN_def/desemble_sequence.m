function out = desemble_sequence(seq_aligned, ref_aligned)

    chop_sites = CARLIN_def.partition_sequence(seq_aligned, ref_aligned);
    
    seq_segments = arrayfun(@(i) seq_aligned(chop_sites(i,1):chop_sites(i,2)), 1:size(chop_sites,1), 'UniformOutput', false);
    ref_segments = arrayfun(@(i) ref_aligned(chop_sites(i,1):chop_sites(i,2)), 1:size(chop_sites,1), 'UniformOutput', false);

    out = AlignedSEQ(seq_segments, ref_segments);
    
    assert(strcmp(seq_aligned, out.get_seq()));
    assert(strcmp(ref_aligned, out.get_ref()));

end
