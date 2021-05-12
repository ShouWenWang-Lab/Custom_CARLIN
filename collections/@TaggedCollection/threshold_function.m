function thresholds = threshold_function(freqs, max_elem, Nreads, p, read_floor, read_override, type)

    if (strcmp(type, 'UMI'))
        field_name = 'max_molecules';
    else
        field_name = 'max_cells';
    end
    
    ind_99_pctl = max(round(size(freqs,1)/100),1);
    thresholds.one_tenth_99_pctl = ceil(freqs(ind_99_pctl)/10);
    thresholds.(field_name) = freqs(min(max_elem, size(freqs,1)))+1;    
    thresholds.equal_partition = ceil(Nreads/max_elem);    
    thresholds.err_floor = ceil(freqs(1)*p*(1-p)^9);
    thresholds.read_floor = read_floor;
    if (isnan(read_override))
        thresholds.override = read_override;
        thresholds.chosen = nanmax(structfun(@(x) x, thresholds));
    else
        thresholds.override = read_override;
        thresholds.chosen = thresholds.override;
    end
    fprintf('...(%d,%d,%d,%d,%d,%d,%d) reads for (99th pctl/10,%s,equal,seq_err,min,override,final)\n', ...
        thresholds.one_tenth_99_pctl, thresholds.(field_name), thresholds.equal_partition, ...
        thresholds.err_floor, thresholds.read_floor, thresholds.override, thresholds.chosen, field_name);
    
end