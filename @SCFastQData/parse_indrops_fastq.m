function [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = parse_indrops_fastq(fastq_file)

    if (ischar(fastq_file))
        fastq_file = cellstr(fastq_file);
    end
        
    assert(size(fastq_file,2) == 1, 'FASTQ files should be supplied as a semi-colon separated list.');
    
    Nfastqs = size(fastq_file,1);
    Nreads = cell(Nfastqs,1);
    
    CB = cell(Nfastqs,1);
    UMI = cell(Nfastqs,1);
    QC  = cell(Nfastqs,1);
    SEQ = cell(Nfastqs,1);
    
    for i = 1:Nfastqs

        assert(exist(fastq_file{i}, 'file') == 2, 'Missing FASTQ file');

        fprintf('Parsing FASTQ file: %s\n', fastq_file{i});
        
        [fastq_file{i}, ext] = FastQData.maybe_unzip(fastq_file{i});
        
        Nreads{i} = fastqinfo(fastq_file{i});
        Nreads{i} = Nreads{i}.NumberOfEntries;
        fprintf('Reads in FASTQ: %d\n', Nreads{i});

        H = cell(Nreads{i},1);    
        [H(:,1), ~, ~] = fastqread(fastq_file{i});    
        [CB{i}, UMI{i}, QC{i}] = SCFastQData.parse_indrops_provenance(H); 
        clear H;

        SEQ{i} = cell(Nreads{i}, 1);
        [~, SEQ{i}(:,1), ~] = fastqread(fastq_file{i});
        
        FastQData.maybe_clear_unzipped(fastq_file{i}, ext);
        
    end
    
    Nreads = sum(cell2mat(Nreads));
    
    if (Nfastqs > 1)
        fprintf('Total reads in FASTQs: %d\n', Nreads);
    end
    
    QC = vertcat(QC{:});
    CB = vertcat(CB{:});
    UMI = vertcat(UMI{:});
    SEQ = vertcat(SEQ{:});
    
    [CB, ~, read_CB]   = unique_by_freq(CB);
    [UMI, ~, read_UMI] = unique_by_freq(UMI);
    [SEQ, ~, read_SEQ] = unique_by_freq(SEQ);
    
end