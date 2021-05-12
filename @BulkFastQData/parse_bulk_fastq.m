function [SEQ, read_SEQ, QC, Nreads] = parse_bulk_fastq(fastq_file, cfg)

    if (ischar(fastq_file))
        fastq_file = cellstr(fastq_file);
    end
        
    assert(size(fastq_file,2) == 1, 'FASTQ files should be supplied as a semi-colon separated list.');
    
    Nfastqs = size(fastq_file,1);
    Nreads = cell(Nfastqs,1);
    SEQ = cell(Nfastqs,1);
    QC  = cell(Nfastqs,1);
    
    for i = 1:Nfastqs

        assert(exist(fastq_file{i}, 'file') == 2, sprintf('Missing FASTQ file: %s', fastq_file{i}));
        
        fprintf('Parsing FASTQ file: %s\n', fastq_file{i});

        [fastq_file{i}, ext] = FastQData.maybe_unzip(fastq_file{i});
        
        Nreads{i} = fastqinfo(fastq_file{i});
        Nreads{i} = Nreads{i}.NumberOfEntries;
        fprintf('Reads in FASTQ: %d\n', Nreads{i});

        SEQ{i} = cell(Nreads{i},1);
        QC{i} = cell(Nreads{i},1);

        [~, SEQ{i}(:,1), QC{i}(:,1)] = fastqread(fastq_file{i});

        FastQData.maybe_clear_unzipped(fastq_file{i}, ext);
    end
    
    Nreads = sum(cell2mat(Nreads));
    
    if (Nfastqs > 1)
        fprintf('Total reads in FASTQs: %d\n', Nreads);
    end
    
    SEQ = vertcat(SEQ{:});
    QC  = vertcat(QC{:});
    
    [SEQ, ~, read_SEQ] = unique_by_freq(SEQ);
    [SEQ, QC] = FastQData.orient_reads(cfg, SEQ, QC);
    
end