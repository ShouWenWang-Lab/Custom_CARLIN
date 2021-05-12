function [CB, UMI, QC] = parse_10x_provenance(CB, QC, cfg)
    
    QC = cellfun(@(x) x(1:min(length(x),cfg.CB.length+cfg.UMI.length)), QC, 'un', false);
    CB = cellfun(@(x) x(1:min(length(x),cfg.CB.length+cfg.UMI.length)), CB, 'un', false);
    [CB, UMI] = cellfun(@(x) deal(x(1:min(length(x),cfg.CB.length)), x(min(length(x),cfg.CB.length)+1:end)), CB, 'un', false);
    
end