classdef (Abstract) TaggedCollection < handle
    
    methods (Access = public)
        [out, collapse_map] = denoise(obj, vargin);
        thresholds = compute_thresholds(obj, params, FQ, varargin);
    end
    
    methods (Static)    
        tag_map = directional_adjacency_top_down_denoiser(tags, weights, exclude);
        thresholds = threshold_function(freqs, max_elem, Nreads, p, read_floor, read_override, type)
    end
end