function masks = filter_sc_CBs_and_UMIs(cfg, CB, read_CB, UMI, read_UMI, QC)

    fprintf('Analyzing SC CBs and UMIs\n');

    assert(strcmp(cfg.type, 'SC'), 'Invalid CFG passed function');
    assert(size(read_CB,1) == size(QC,1), 'Dimension mismatch QC and CB');
    assert(size(read_UMI,1) == size(QC,1), 'Dimension mismatch QC and UMI');

    masks.QC_length_match = cellfun(@length, CB(read_CB)) + cellfun(@length, UMI(read_UMI)) == cellfun(@length, QC);

    if(strcmp(cfg.SC.Platform, 'InDrops'))        
        if (cfg.SC.Version == 3)
            masks.CB_correct_length = cellfun(@(x) length(x) == cfg.CB(1).length+cfg.CB(2).length, CB);
        else
            masks.CB_correct_length = cellfun(@(x) length(x) >= cfg.CB(1).length(1)+cfg.CB(2).length && ...
                                                   length(x) <= cfg.CB(1).length(2)+cfg.CB(2).length, CB);
        end
        masks.UMI_correct_length = cellfun(@(x) length(x) == cfg.UMI.length, UMI);        
    elseif(strcmp(cfg.SC.Platform, '10x'))
        masks.CB_correct_length = cellfun(@(x) length(x) == cfg.CB.length, CB);
        masks.UMI_correct_length = cellfun(@(x) length(x) == cfg.UMI.length, UMI);        
    end

    masks.CB_correct_length = masks.CB_correct_length(read_CB);
    masks.UMI_correct_length = masks.UMI_correct_length(read_UMI);    

    masks.CB_no_N = ~cellfun(@(x) any(x=='N'), CB);
    masks.CB_no_N = masks.CB_no_N(read_CB);
    
    masks.UMI_no_N = ~cellfun(@(x) any(x=='N'), UMI);
    masks.UMI_no_N = masks.UMI_no_N(read_UMI);
    
    masks.good_CB_UMI_QC = cellfun(@(x) all(double(x)-33 >= 20), QC);

    masks.valid_provenance_structure = ( masks.QC_length_match    & masks.CB_correct_length & masks.UMI_correct_length & ...
                                         masks.good_CB_UMI_QC     & masks.CB_no_N           & masks.UMI_no_N );

    masks.QC_length_match            = uint32(find(masks.QC_length_match));
    masks.CB_correct_length          = uint32(find(masks.CB_correct_length));
    masks.UMI_correct_length         = uint32(find(masks.UMI_correct_length));
    masks.good_CB_UMI_QC             = uint32(find(masks.good_CB_UMI_QC));
    masks.CB_no_N                    = uint32(find(masks.CB_no_N));
    masks.UMI_no_N                   = uint32(find(masks.UMI_no_N));
    masks.valid_provenance_structure = uint32(find(masks.valid_provenance_structure));

    N = length(read_CB);

    fprintf('From %d reads, found good (L_QC,L_CB,L_UMI,QC,no N CB,no N UMI,all) (%d,%d,%d,%d,%d,%d,%d) times\n', ...
        N, length(masks.QC_length_match), length(masks.CB_correct_length), length(masks.UMI_correct_length), ...
        length(masks.good_CB_UMI_QC), length(masks.CB_no_N), length(masks.UMI_no_N), length(masks.valid_provenance_structure));

end