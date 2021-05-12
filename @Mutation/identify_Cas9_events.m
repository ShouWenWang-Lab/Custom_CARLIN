function event_list = identify_Cas9_events(aligned)

    assert(isa(aligned, 'AlignedSEQ'));

    [mut_list, ~, seq, ref_seq, ref_mask] = Mutation.identify_sequence_events(aligned);
    ref = CARLIN_def.getInstance;
    
    if (~isempty(mut_list))        
        compound = zeros(0,2);
        compound_start = 1;
        compound_end = 1;
        for i = 2:size(mut_list,1)        
            start_loc_bp = mut_list(i).loc_start;
            end_loc_bp   = mut_list(compound_end).loc_end;            
            so_far_site = CARLIN_def.locate(ref, end_loc_bp , ref.bounds.ordered);
            new_site    = CARLIN_def.locate(ref, start_loc_bp, ref.bounds.ordered);        
            if ((start_loc_bp == end_loc_bp+1) || (new_site.abs == so_far_site.abs) || ...
                (new_site.abs==so_far_site.abs+1 && strcmp(so_far_site.type, 'cutsites')) || ...
                start_loc_bp < end_loc_bp)
                compound_end = i;
            else
                compound = [compound; [compound_start compound_end]];
                compound_start = i;
                compound_end = i;
            end
        end       
        compound = [compound; compound_start compound_end];
        event_list = cell(size(compound,1),1);        
        for i = 1:size(event_list,1)
            if (compound(i,1) == compound(i,2))
                event_list{i} = mut_list(compound(i,1));
                % Trivial case of a single insertion, but sequence events
                % don't treat insertions as compound, so patch that up
                % here.
                if (event_list{i}.type == 'I')
                    if (event_list{i}.seq_new(1) ~= event_list{i}.seq_old(1) && ...
                        event_list{i}.seq_new(end) ~= event_list{i}.seq_old(end))
                        event_list{i} = Mutation('C', event_list{i}.loc_start, event_list{i}.loc_end, event_list{i}.seq_new, event_list{i}.seq_old);
                    end
                end
            else
                ref_start_pos = mut_list(compound(i,1)).loc_start;
                temp_list = mut_list(compound(i,1):compound(i,2));
                ref_end_pos   = max([temp_list.loc_end]);
                seq_start_pos = ref_mask(ref_start_pos);
                seq_end_pos   = ref_mask(ref_end_pos);
                if (mut_list(compound(i,1)).type == 'I' && mut_list(compound(i,1)).seq_old(1)=='-')
                    if (ref_start_pos > 1)
                        seq_start_pos = ref_mask(ref_start_pos-1)+1;
                    else
                        seq_start_pos = 1;
                    end
                end
                if ((mut_list(compound(i,2)).type == 'I' && mut_list(compound(i,2)).seq_old(end)=='-') && ...
                    mut_list(compound(i,2)).loc_end >= mut_list(compound(i,1)).loc_end)
                    if (ref_end_pos < length(ref_mask))
                        seq_end_pos = ref_mask(ref_end_pos+1)-1;
                    else
                        seq_end_pos = length(seq);
                    end
                end
                if (mut_list(compound(i,1)).type == 'M')
                    seq_start_pos = max(seq_start_pos-1,1);
                    ref_start_pos = max(ref_start_pos-1,1);
                end
                if (mut_list(compound(i,2)).type == 'M')
                    seq_end_pos = min(seq_end_pos+1, length(seq));
                    ref_end_pos = min(ref_end_pos+1,length(ref_mask));
                end 
                compound_ref_seq = ref_seq(ref_mask(ref_start_pos):ref_mask(ref_end_pos));
                compound_in_seq = seq(seq_start_pos:seq_end_pos);
                event_list{i} = Mutation.make_compound_mutation(ref_start_pos, ref_end_pos, compound_in_seq, compound_ref_seq);
            end
        end
        
        event_list = vertcat(event_list{:});
        
        mutated_allele = Mutation.apply_mutations(event_list);

        assert(isequal(degap(mutated_allele.get_seq), degap(seq)), '%s\n%s\n', degap(mutated_allele.get_seq), degap(seq));
        assert(isequal(degap(mutated_allele.get_ref), degap(ref_seq)));
    else
        event_list = [];
    end    
end