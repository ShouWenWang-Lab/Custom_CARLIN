 function masks = filter_bulk_UMIs(UMI, read_UMI, QC, L,QC_threshold)
   % SW: I updated the QC_threshold choice. If this is not provided, reverted to the original method
    if nargin<5
      QC_threshold = 30; % if not specified, use the default
      orig_method=true;
    else
      orig_method=false;
    end 
   
            
    fprintf('Analyzing Bulk UMIs\n');

    assert(L > 0, 'Expected UMI length (%d) must be positive', L);
    assert(size(read_UMI,1) == size(QC,1), 'Mismatch in dimension of UMI and QC');
    assert(isequal(cellfun(@length, UMI(read_UMI)), cellfun(@length, QC)), 'QC and UMI lengths do not match');
    assert(all(cellfun(@length, UMI) <= L), 'Supplied UMIs should not be greater than supplied length');

    masks.UMI_correct_length = cellfun(@length, UMI) == L;
    masks.UMI_correct_length = masks.UMI_correct_length(read_UMI);

    if orig_method
        masks.UMI_good_QC = cellfun(@(x) all(double(x)-33 >= QC_threshold), QC); %QC score control 
    else
        masks.UMI_good_QC = cellfun(@(x) mean(double(x)-33) >= QC_threshold, QC); %QC score control 
    end

    masks.UMI_no_N = ~cellfun(@(x) any(x=='N'), UMI);
    masks.UMI_no_N = masks.UMI_no_N(read_UMI);

    masks.valid_provenance_structure = masks.UMI_correct_length & masks.UMI_good_QC & masks.UMI_no_N;

    masks.UMI_correct_length = uint32(find(masks.UMI_correct_length));
    masks.UMI_good_QC        = uint32(find(masks.UMI_good_QC));
    masks.UMI_no_N           = uint32(find(masks.UMI_no_N));
    masks.valid_provenance_structure = uint32(find(masks.valid_provenance_structure));

    fprintf('From %d reads, found (length,QC,no N,all) UMIs (%d,%d,%d,%d) times\n', ...
        length(read_UMI), length(masks.UMI_correct_length), length(masks.UMI_good_QC), ...
        length(masks.UMI_no_N), length(masks.valid_provenance_structure));
end