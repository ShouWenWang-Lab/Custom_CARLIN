function event = classify_motif_event(seq, ref)

    assert(length(seq) == length(ref));
    assert(ischar(seq) && ischar(ref));
    
    if (strcmp(seq, ref))
        event = 'N';
    elseif (all(seq=='-'))
        event = 'E';
    elseif (all(ref~='-'))
        if (any(seq == '-'))
            event = 'D';
        else
            event = 'M';
        end
    else
        event = 'I';
    end
    
end