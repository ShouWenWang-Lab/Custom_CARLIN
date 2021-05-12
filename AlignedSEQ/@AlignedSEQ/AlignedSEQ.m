classdef AlignedSEQ
    
    properties (SetAccess = immutable)
        motifs
    end
    
    methods (Access = public)
        
        function obj = AlignedSEQ(seq_segments, ref_segments)
            assert(length(seq_segments) == length(ref_segments));
            assert(length(seq_segments) == CARLIN_def.getInstance.N.motifs);
            obj.motifs = cellfun(@(s,r) AlignedSEQMotif(s, r), seq_segments, ref_segments, 'un', false);
            obj.motifs = vertcat(obj.motifs{:});
        end
        
        function seq = get_seq(obj)
            seq = [obj.motifs.seq];
        end
        
        function ref = get_ref(obj)
            ref = [obj.motifs.ref];
        end
        
        function event = get_event_structure(obj)
            event = [obj.motifs.event];
        end
        
    end
    
    methods (Static)
        out = sanitize_conserved_regions(in);
        out = sanitize_prefix_postfix(in);
    end
end