function [CB, UMI, QC] = parse_indrops_provenance(H)

    H = cellfun(@(x) x(1:find(x==':', 1, 'last')-1), H, 'un', false);
    
    cb_halves_split_loc = regexp(H, '-', 'once');
    cb_umi_split_loc = regexp(H, ':');
    qc_split_loc = cellfun(@(x) x(2), cb_umi_split_loc, 'un', false);
    cb_umi_split_loc = cellfun(@(x) x(1), cb_umi_split_loc, 'un', false);

    [CB, UMI, QC] = cellfun(@(x,c,u,q) deal([x(1:c-1) x(c+1:u-1)], x(u+1:q-1), [x(q+1:q+c-1) x(q+c+1:q+u-1) x(q+u+1:q+q-1)]), ...
        H, cb_halves_split_loc, cb_umi_split_loc, qc_split_loc, 'un', false);

end