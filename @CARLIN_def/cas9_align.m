function [sc, al] = cas9_align(seq)

    assert(~isempty(seq), 'Sequence to align is empty');
    assert(~any(isgap(seq)), 'Input sequence should not have any insertions');
    
    ref           = CARLIN_def.getInstance;
    refseq        = ref.seq.CARLIN;
    open_penalty  = ref.open_penalty;
    close_penalty = ref.close_penalty;
    match         = ref.sub_score;
    
    [sc, al] = CARLIN_def.cas9_align_mex([0 nt2int(seq)], [0 nt2int(refseq)], open_penalty, [0 close_penalty], padarray(match, [1 1], 'pre'));        
    al(al==0)=16;
    al = int2nt(al);
    al_seq = al(1,:);
    al_ref = al(2,:);
   
    assert(isequal(degap(al_seq), seq), '%s\n%s', degap(al_seq), seq);
    assert(isequal(degap(al_ref), refseq), '%s\n%s', degap(al_ref), refseq);
    assert(length(al_seq) == length(al_ref));
    assert(~any(isgap(al_seq) & isgap(al_ref)));
    
    al = CARLIN_def.desemble_sequence(al_seq, al_ref);
    
end