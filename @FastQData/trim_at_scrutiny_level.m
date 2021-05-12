function [valid_mask, ind] = trim_at_scrutiny_level(level, SEQ, trim_SEQ, which_end, offset, thresh)

    switch level
        
        case 'exact' % only cares the accuracy of start or end position
            
            ind = strfind(SEQ, trim_SEQ); %the starting indices of any occurrences of trim_SEQ in SEQ
            valid_mask = ~cellfun(@isempty, ind);
            if (strcmp(which_end, 'head'))
                valid_mask(valid_mask) = cellfun(@(x) x(1)==offset+1, ind(valid_mask));
                ind(valid_mask) = cellfun(@(x) x(1)+length(trim_SEQ), ind(valid_mask), 'un', false');
            else
                valid_mask(valid_mask) = cellfun(@(x,y) x(end)==length(y)-length(trim_SEQ)-offset+1, ...
                                                             ind(valid_mask), SEQ(valid_mask));
                ind(valid_mask) = cellfun(@(x) x(end)-1, ind(valid_mask), 'un', false');
            end
            
        case 'malformed' % cares about whether the alignment is right or not.
            
            N = size(SEQ,1);
            sc = zeros(N,1);
            al = cell(N,1);
            parfor i = 1:N
                if (~isempty(SEQ{i}))
                    [sc(i), al{i}] = nwalign(trim_SEQ, SEQ{i}, 'Alphabet', 'NT', 'Glocal', true);
                end
            end
            valid_mask = sc >= thresh;
            ind = cell(1,N);
            if (strcmp(which_end, 'head'))
                ind(valid_mask) = cellfun(@(x) find(x(1,:)~='-', 1, 'last'), al(valid_mask), 'UniformOutput', false);
                ind(valid_mask) = cellfun(@(x,s) sum(x(3,1:s)~='-')+1, al(valid_mask), ind(valid_mask)', 'UniformOutput', false);
            else
                ind(valid_mask) = cellfun(@(x) find(x(1,:)~='-', 1, 'first'), al(valid_mask), 'UniformOutput', false);
                ind(valid_mask) = cellfun(@(x,e) sum(x(3,1:e-1)~='-'), al(valid_mask), ind(valid_mask)', 'UniformOutput', false);
            end
            ind = ind';
            
        case 'misplaced'
            ind = strfind(SEQ, trim_SEQ);
            valid_mask = ~cellfun(@isempty, ind);
            if (strcmp(which_end, 'head'))                
                ind(valid_mask) = cellfun(@(x) x(1)+length(trim_SEQ), ind(valid_mask), 'un', false');
            else
                ind(valid_mask) = cellfun(@(x) x(end)-1, ind(valid_mask), 'un', false');
            end
            
        case 'ignore'
            valid_mask = true(size(SEQ,1),1);
            if (strcmp(which_end, 'head'))
                ind = repmat({1}, size(SEQ,1), 1);
            else
                ind = cellfun(@length, SEQ, 'UniformOutput', false);
            end
            
        otherwise
            error('Unrecognized trim level: %s', level);
    end
    ind(~valid_mask) = {0};
    ind = uint16(vertcat(ind{:}));
end
    