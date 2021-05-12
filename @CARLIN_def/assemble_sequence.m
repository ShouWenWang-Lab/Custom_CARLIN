function out = assemble_sequence(segments)

    segments = order_named_motifs(segments);
    out = AlignedSEQ(segments(:,1), segments(:,2));
    
end
