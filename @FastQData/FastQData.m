classdef (Abstract=true) FastQData < handle
    
    properties (SetAccess = immutable, GetAccess = public)
        source_files
        Nreads
        SEQ_raw
        read_SEQ_raw
        SEQ_trimmed
        read_SEQ_trimmed
        SEQ_valid
        read_SEQ_valid
        QC
        UMI
        read_UMI
        masks
        trim_loc
    end
    
        
    methods (Static)
        
        % Zip/unzip helper
        [fastq_file, ext] = maybe_unzip(fastq_file);
        maybe_clear_unzipped(fastq_file, ext);
        
        % Trim CARLIN from read
        [SEQ, QC] = orient_reads(cfg, SEQ, QC)
        [valid_mask, ind] = trim_at_scrutiny_level(level, SEQ, trim_SEQ, which_end, offset, thresh);
        [SEQ_trimmed, read_SEQ_trimmed, masks, trim_loc] = extract_CARLIN_from_sequences(SEQ_raw, read_SEQ_raw, cfg);
        
    end
    
    methods (Access = public)
        
        function obj = FastQData(source_files, Nreads, SEQ_raw, read_SEQ_raw, SEQ_trimmed, read_SEQ_trimmed, ...
                                 UMI, read_UMI, QC, masks, trim_loc)
                             
            obj.source_files = source_files;
            obj.Nreads = uint32(Nreads);
            obj.SEQ_raw = cellfun(@nt2int, SEQ_raw, 'un', false);
            obj.read_SEQ_raw = uint32(read_SEQ_raw);
            obj.SEQ_trimmed = cellfun(@nt2int, SEQ_trimmed, 'un', false);
            obj.read_SEQ_trimmed = uint32(read_SEQ_trimmed);
            
            [r, ~, v] = find(accumarray(obj.read_SEQ_trimmed(masks.valid_lines),1));
            [~, idx] = sort(v, 'descend');
            r = r(idx);
            obj.SEQ_valid = obj.SEQ_trimmed(r);
            obj.read_SEQ_valid = uint32(zeros(size(obj.read_SEQ_trimmed)));
            [~, obj.read_SEQ_valid(masks.valid_lines)] = ismember(obj.read_SEQ_trimmed(masks.valid_lines), r);
            
            obj.UMI = UMI;
            obj.read_UMI = uint32(read_UMI);
            obj.QC = cellfun(@uint8, QC, 'un', false);
            obj.masks = masks;
            obj.trim_loc = trim_loc;
            
        end
        
        function [SEQ_ind, SEQ_weight] = get_SEQ_ind_by_UMI(obj, UMI, filter)
            
            if (size(UMI,1) > 1)                
                assert(iscell(UMI));
            else
                UMI = cell(UMI);
            end
            
            if (nargin == 3)
                filter = intersect(obj.masks.valid_lines, filter);
            else
                filter = obj.masks.valid_lines;
            end
            
            [is, which_UMI] = ismember(UMI, obj.UMI);            
            if (any(is))
                filter = filter(ismember(obj.read_UMI(filter), which_UMI(is)));
                [SEQ_ind, ~, SEQ_weight] = unique([obj.read_UMI(filter) obj.read_SEQ_valid(filter)], 'rows');
                SEQ_weight = accumarray(SEQ_weight,1);
                [SEQ_ind, SEQ_weight] = ...
                    arrayfun(@(i) deal(SEQ_ind(SEQ_ind(:,1)==i,2), SEQ_weight(SEQ_ind(:,1)==i)), which_UMI(is), 'un', false);
                SEQ_ind = cellfun(@(x,y) sortrows([x, y], 2, 'descend'), SEQ_ind, SEQ_weight, 'un', false);
                [SEQ_ind, SEQ_weight] = cellfun(@(x) deal(x(:,1), x(:,2)), SEQ_ind, 'un', false);
                SEQ_ind(is) = SEQ_ind;
                SEQ_weight(is) = SEQ_weight;
            end
            
            SEQ_ind(~is) = {[]};
            SEQ_weight(~is) = {[]};
            
        end
        
        function UMIs = get_UMIs(obj)
            UMIs = obj.UMI(unique_by_freq(obj.read_UMI(obj.masks.valid_lines)));
        end
        
        function SEQs = get_SEQs(obj)
            assert(~any(cellfun(@isempty, obj.SEQ_valid)));
            SEQs = cellfun(@int2nt, obj.SEQ_valid, 'un', false);
        end
    end

end
