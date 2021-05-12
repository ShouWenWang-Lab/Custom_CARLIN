function RGB = get_sequence_coloring(aligned_seqs, resolution)

    assert(iscell(aligned_seqs) && isa(aligned_seqs{1}, 'AlignedSEQ'));
    if (nargin < 2)
        resolution = 'bp';
    end
    
    ref = CARLIN_def.getInstance;
    
    switch(resolution)
        case 'bp'
            event = cellfun(@(x) Mutation.classify_bp_event(x), aligned_seqs, 'Un', false);
            event = vertcat(event{:});
        case 'motif'
            event = cellfun(@(x) AlignedSEQMotif.classify_motif_event(x.get_seq(), x.get_ref()), aligned_seqs, 'Un', false);
            event = vertcat(event{:});
            event = repelem(event, 1, diff(ref.bounds.ordered,[],2)+1);
        otherwise
            error("Unrecognized sequence coloring resolution");
    end
    
    RGB = arrayfun(@(x) ref.color(x), event, 'un', false);
    RGB = reshape(vertcat(RGB{:}), [], ref.width.CARLIN, 3);
    
end
