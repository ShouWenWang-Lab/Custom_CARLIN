function aligned_seqs = sanitize_conserved_regions(aligned_seqs)

    is_obj = false;
    if (isa(aligned_seqs, 'AlignedSEQ'))
        is_obj = true;
        aligned_seqs = {aligned_seqs};
    end

    assert(iscell(aligned_seqs) && isa(aligned_seqs{1}, 'AlignedSEQ'));
    
    L_start = cellfun(@(x) length(x.get_seq()), aligned_seqs);
    
    seq_events = cellfun(@(x) x.get_event_structure(), aligned_seqs, 'un', false);
    seq_events = vertcat(seq_events{:});
    ref = CARLIN_def.getInstance;

    motif_candidate = true(1, ref.N.motifs);
    motif_candidate(ref.motifs.cutsites) = false;    
    
    motifs_to_clean = (seq_events == 'M') & motif_candidate;
    seqs_to_clean = any(motifs_to_clean,2);
    
    if (any(seqs_to_clean))
    
        [seq_motifs, ref_motifs] = cellfun(@(x) deal({x.motifs.seq}, {x.motifs.ref}), aligned_seqs(seqs_to_clean), 'un', false);

        motifs_to_clean = find(motifs_to_clean(seqs_to_clean,:));
        seq_motifs = vertcat(seq_motifs{:});
        ref_motifs = vertcat(ref_motifs{:});

        assert(~any(cellfun(@(x) any(x=='-'), ref_motifs(motifs_to_clean))));
        assert(~any(cellfun(@(x) any(x=='-'), seq_motifs(motifs_to_clean))));
        
        seq_motifs(motifs_to_clean) = ref_motifs(motifs_to_clean);

        seq_motifs = arrayfun(@(i) seq_motifs(i,:), [1:size(seq_motifs,1)]', 'un', false);
        ref_motifs = arrayfun(@(i) ref_motifs(i,:), [1:size(ref_motifs,1)]', 'un', false);   

        aligned_seqs(seqs_to_clean) = cellfun(@(s,r) AlignedSEQ(s,r), seq_motifs, ref_motifs, 'un', false);    

        L_end = cellfun(@(x) length(x.get_seq()), aligned_seqs);
        assert(all(L_start==L_end));
    end
    
    if (is_obj)
        aligned_seqs = aligned_seqs{1};
        assert(isa(aligned_seqs, 'AlignedSEQ'));
    end    
end