function [uniq_val, backtrack, uniq_ind] = unique_by_freq(val, weights)

    if (nargin < 2)
        weights = 1;
    end
    if (iscell(val))
        [uniq_val, backtrack, uniq_ind] = unique(val);
    else
        [uniq_val, backtrack, uniq_ind] = unique(val, 'rows');
    end
    backtrack = uint32(backtrack);
    uniq_ind = uint32(uniq_ind);
    [~, idx] = sort(accumarray(uniq_ind,weights), 'descend');
    uniq_val = uniq_val(idx,:);
    backtrack = backtrack(idx,:);
    [~, uniq_ind] = ismember(uniq_ind, idx);
end