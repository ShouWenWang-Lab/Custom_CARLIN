function A = merge_structs(A, B)

    fn_A = fieldnames(A);
    fn_B = fieldnames(B);
    assert(~any(ismember(fn_A, fn_B)), 'Cannot merge structures with overlapping field names');
    for i = 1:length(fn_B)
        A.(fn_B{i}) = B.(fn_B{i});
    end
end