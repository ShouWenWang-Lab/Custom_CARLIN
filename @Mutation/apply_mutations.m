function mutated_allele = apply_mutations(event_list)
    
    ref = CARLIN_def.getInstance.seq.CARLIN;
    seq = ref;
    offset = 0;
    
    if (~isempty(event_list))
        assert(isa(event_list, 'Mutation'));
        for i = 1:size(event_list,1)
            seq = [seq(1:event_list(i).loc_start-1+offset) event_list(i).seq_new seq(event_list(i).loc_end+1+offset:end)];            
            ref = [ref(1:event_list(i).loc_start-1+offset) event_list(i).seq_old ref(event_list(i).loc_end+1+offset:end)];            
            if (event_list(i).type=='C' || event_list(i).type=='I')                
                offset = offset+length(event_list(i).seq_new) - (event_list(i).loc_end-event_list(i).loc_start+1);
            end            
        end 
    end
    mutated_allele = CARLIN_def.desemble_sequence(seq, ref);
    
end