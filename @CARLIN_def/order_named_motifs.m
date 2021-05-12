function ordered = order_named_motifs(ref, unordered)
    if (iscell(unordered.prefix))        
        ordered = cell(ref.N.motifs, size(unordered.prefix,2));        
    else
        ordered = zeros(ref.N.motifs, size(unordered.prefix,2));
    end
    ordered(ref.motifs.prefix,:)   = unordered.prefix;
    ordered(ref.motifs.consites,:) = unordered.consites;
    ordered(ref.motifs.cutsites,:) = unordered.cutsites;
    ordered(ref.motifs.pams,:)     = unordered.pams;
    ordered(ref.motifs.postfix,:)  = unordered.postfix;
end