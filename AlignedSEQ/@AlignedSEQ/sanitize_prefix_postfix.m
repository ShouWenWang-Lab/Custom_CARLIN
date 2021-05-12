function aligned_seqs = sanitize_prefix_postfix(aligned_seqs)

    is_obj = false;
    
    if (isa(aligned_seqs, 'AlignedSEQ'))
        is_obj = true;
        aligned_seqs = {aligned_seqs};
    end
    assert(iscell(aligned_seqs) && isa(aligned_seqs{1}, 'AlignedSEQ'));
    
    seq_events = cellfun(@(x) x.get_event_structure(), aligned_seqs, 'un', false);
    seq_events = vertcat(seq_events{:});
    
    head_insert = find(seq_events(:,1) == 'I');
    
    if (~isempty(head_insert))
        [seq_motifs, ref_motifs] = cellfun(@(x) deal({x.motifs.seq}, {x.motifs.ref}), aligned_seqs(head_insert), 'un', false);
        seq_motifs = vertcat(seq_motifs{:});
        ref_motifs = vertcat(ref_motifs{:});
        head_ind = cellfun(@(x) find(x~='-', 1, 'first'), ref_motifs(:,1));    
        trim_head = (head_ind ~= 1);
        if (any(trim_head))
            head_insert = head_insert(trim_head);
            head_ind = num2cell(head_ind(trim_head));
            seq_motifs = seq_motifs(trim_head,:);
            ref_motifs = ref_motifs(trim_head,:);

            [seq_motifs(:,1), ref_motifs(:,1)] = cellfun(@(s,r,h) deal(s(h:end), r(h:end)), ...
                                                         seq_motifs(:,1), ref_motifs(:,1), head_ind, 'un', false);
                                                     
            seq_motifs = arrayfun(@(i) seq_motifs(i,:), [1:size(seq_motifs,1)]', 'un', false);
            ref_motifs = arrayfun(@(i) ref_motifs(i,:), [1:size(ref_motifs,1)]', 'un', false);
            aligned_seqs(head_insert) = cellfun(@(s,r) AlignedSEQ(s,r), seq_motifs, ref_motifs, 'un', false);
        end
    end
    
    tail_insert = find(seq_events(:,end) == 'I');
    
    if (~isempty(tail_insert))
        [seq_motifs, ref_motifs] = cellfun(@(x) deal({x.motifs.seq}, {x.motifs.ref}), aligned_seqs(tail_insert), 'un', false);
        seq_motifs = vertcat(seq_motifs{:});
        ref_motifs = vertcat(ref_motifs{:});
        tail_ind = cellfun(@(x) find(x~='-', 1, 'last'), ref_motifs(:,end));    
        trim_tail = (tail_ind ~= cellfun(@length, ref_motifs(:,end)));
        if (any(trim_tail))
            tail_insert = tail_insert(trim_tail);
            tail_ind = num2cell(tail_ind(trim_tail));
            seq_motifs = seq_motifs(trim_tail,:);
            ref_motifs = ref_motifs(trim_tail,:);

            [seq_motifs(:,end), ref_motifs(:,end)] = cellfun(@(s,r,t) deal(s(1:t), r(1:t)), ...
                                                         seq_motifs(:,end), ref_motifs(:,end), tail_ind, 'un', false);
                                                     
            seq_motifs = arrayfun(@(i) seq_motifs(i,:), [1:size(seq_motifs,1)]', 'un', false);
            ref_motifs = arrayfun(@(i) ref_motifs(i,:), [1:size(ref_motifs,1)]', 'un', false);                                                     
            aligned_seqs(tail_insert) = cellfun(@(s,r) AlignedSEQ(s,r), seq_motifs, ref_motifs, 'un', false);
        end
    end
    
    if (is_obj)
        aligned_seqs = aligned_seqs{1};
        assert(isa(aligned_seqs, 'AlignedSEQ'));
    end
    
end