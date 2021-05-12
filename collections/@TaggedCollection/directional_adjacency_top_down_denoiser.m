function tag_map = directional_adjacency_top_down_denoiser(tags, weights, exclude)

    N = size(tags,1);
    assert(size(weights,1)==N);
    if (nargin < 3)
        exclude = false(N,1);
    else
        assert(size(exclude,1)==N);
    end

    tag_map = zeros(N,1);
    L = cellfun(@length, tags);
    [~, idx] = sort(weights, 'descend');

    for i = idx'
        if (tag_map(i) == 0)
            tag_map(i) = i;
        end
        children = find(weights(i)>=2*weights-1 & tag_map==0 & L==L(i) & ~exclude);
        if (~isempty(children))
            children = children(sum(vertcat(tags{children})~=tags{i},2)==1);
            tag_map(children) = tag_map(i);
        end
    end
end