classdef AlignedSEQMotif
    
    properties (SetAccess = private, GetAccess = public)
        seq
        ref
        event
    end
    
    methods (Access = public)
        
        function obj = AlignedSEQMotif(seq, ref)
            assert(length(seq) == length(ref), ...
                'Sequence length not equal to reference length in motif initialization: %d =/= %d', length(seq), length(ref))
            obj.seq = seq;
            obj.ref = ref;
            obj.event = AlignedSEQMotif.classify_motif_event(seq, ref);
        end
        
    end
    
    methods (Static)
        
        event = classify_motif_event(seq, ref);
        
    end
end
        
        
        