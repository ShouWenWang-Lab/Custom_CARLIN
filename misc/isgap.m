function seq = isgap(seq)
    if (ischar(seq))
        seq = (seq == '-');
    elseif (isinteger(seq))
        seq = (seq == 16);
    else
        error('Unknown sequence type'); 
    end
end