function seq = degap(seq)
    seq = seq(~isgap(seq));
end