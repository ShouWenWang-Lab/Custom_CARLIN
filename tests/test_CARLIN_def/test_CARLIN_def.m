classdef (Sealed) test_CARLIN_def < handle        
    %% CARLIN Sequence
    properties (Constant, GetAccess=public) % that this should be public access
        
    
    %% change only below
    
        segments = {'GCCGGCGAGCGCTATGAGCGACT';
            'AGTCGACACGACTCGCGCATACG';
            'AGTCGACTACAGTCGCTACGACG';
            'AGTCGATACGATACGCGCACGCT';
            'AGTCGACTGCACGACAGTCGACG';
            'AGTCGATACGTAGCACGCAGACG';
            'GAGCGAGTCGAGACGCTGACGAT';
            'AGTCGATAGTATGCGTACACGCG';
            'AGTCGCGACTGTACGCACACGCG';
            'AGTCGAGAGCGCGCTCGTCGACT'};
        prefix = 'GC';
        pam = 'ATGG';
        postfix = 'A';
        Primer5 = 'CCTAGCCGGGGATCCTCTAGAGTCGAATGTACAAGTAAAGCGGCC';
        Primer3 = 'TCTAGTTGC';
        SecondarySequence = 'GTCTGCTGTGTGCCT';
        
%         segments = {'TAGTCGACGAGCGCGCTCTCGCTC';
%                     'TCGTCTGCGTGCTACGTATCGACT';
%                     'TCGTCGACTGTCGTGCAGTCGACT';
%                     'TCGCGTGTGCGTACAGTCGCGACT';
%                     'TAGCGTGCGCGTATCGTATCGACT';
%                     'TCGTCGTAGCGACTGTAGTCGACT';
%                     'TCGCGTGTACGCATACTATCGACT';
%                     'TAGTCGCTCATAGCGCTCGCGACT';
%                     'TCGTATGCGCGAGTCGTGTCGACT';
%                     'TATCGTCAGCGTCTCGACTCCGGC'};
%         prefix = 'AGACT';%between Primer5 and first CARLIN segment, part of CARLIN sequnece
%         pam = 'CCA'; 
%         postfix = 'G'; % this is not right. Does it matter?, part of CARLIN sequnece
%         Primer5 = 'GCAACTAGAAGGCACCGACA'; % the very beginning, TAGAAGGCACCGACA, not part of CARLIN SEQ 
%         Primer3 = 'GC'; %not part of CARLIN SEQ, need to reach the very end, not used to quality control (ignored). The last bases tend to be noisy  
%         SecondarySequence = 'ATTCGCGAGGTACCGA'; %after postfix, and before Primer3, not part of CARLIN SEQ, used to ensure the whole sequence is read 
    %% stop change
    
    
    
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

