classdef (Sealed) test_CARLIN_def < handle        
    %% CARLIN Sequence
    properties (Constant, GetAccess=public)
        segments = {'TAGTCGACGAGCGCGCTCTCGCTC';
                    'TCGTCTGCGTGCTACGTATCGACT';
                    'TCGTCGACTGTCGTGCAGTCGACT';
                    'TCGCGTGTGCGTACAGTCGCGACT';
                    'TAGCGTGCGCGTATCGTATCGACT';
                    'TCGTCGTAGCGACTGTAGTCGACT';
                    'TCGCGTGTACGCATACTATCGACT';
                    'TAGTCGCTCATAGCGCTCGCGACT';
                    'TCGTATGCGCGAGTCGTGTCGACT';
                    'TATCGTCAGCGTCTCGACTCCGGC'};
        prefix = 'AGACT';%between Primer5 and first CARLIN segment, part of CARLIN sequnece
        pam = 'CCA'; 
        postfix = 'G'; % this is not right. Does it matter?, part of CARLIN sequnece
        Primer5 = 'GCAACTAGAAGGCACCGACA'; % the very beginning, TAGAAGGCACCGACA, not part of CARLIN SEQ 
        Primer3 = 'GC'; %not part of CARLIN SEQ, need to reach the very end, not used to quality control (ignored). The last bases tend to be noisy  
        SecondarySequence = 'ATTCGCGAGGTACCGA'; %after postfix, and before Primer3, not part of CARLIN SEQ, used to ensure the whole sequence is read 
    end
    
        %% Constant Definitions
    properties (SetAccess = public)
        N;
        motifs;
        seq;
        width;
        bounds;
        bps;
        color;
        alpha;
        match_score;
        open_penalty;
        close_penalty;
        sub_score;
    end 
    
end

