function [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = parse_10x_fastq(fastq_file, cfg)

    assert(size(fastq_file,2)==2, '10x analysis needs comma separated pairs of FASTQ files');
    
    Nfastqs = size(fastq_file,1);
    Nreads = cell(Nfastqs,2);
    
    CB = cell(Nfastqs,1);
    UMI = cell(Nfastqs,1);
    QC  = cell(Nfastqs,1);
    SEQ = cell(Nfastqs,1);
    
    for i = 1:Nfastqs
        
        assert(~strcmp(fastq_file{i,1}, fastq_file{i,2}), ...
            sprintf('Same FASTQ supplied for both CB/UMI and Sequences: %s', fastq_file{i,1}));

        assert(exist(fastq_file{i,1}, 'file') == 2, sprintf('Missing FASTQ file: %s', fastq_file{i,1}));
        assert(exist(fastq_file{i,2}, 'file') == 2, sprintf('Missing FASTQ file: %s', fastq_file{i,2}));
        
        fprintf('Parsing FASTQ files: %s, %s\n', fastq_file{i,1}, fastq_file{i,2});

        [fastq_file(i,:), ext] = FastQData.maybe_unzip(fastq_file(i,:));
        
        Nreads{i,1} = fastqinfo(fastq_file{i,1});
        Nreads{i,1} = Nreads{i,1}.NumberOfEntries;
        Nreads{i,2} = fastqinfo(fastq_file{i,2});
        Nreads{i,2} = Nreads{i,2}.NumberOfEntries;
        assert(Nreads{i,1} == Nreads{i,2}, 'Two 10x FASTQ files have different number of lines');

        fprintf('Reads in FASTQ: %d\n', Nreads{i,1});    

        CB{i} = cell(Nreads{i,1},1);
        QC{i} = cell(Nreads{i,1},1);

        [~, CB{i}(:,1), QC{i}(:,1)] = fastqread(fastq_file{i,1});
        [CB{i}, UMI{i}, QC{i}] = SCFastQData.parse_10x_provenance(CB{i}, QC{i}, cfg);

        SEQ{i} = cell(Nreads{i,1},1);
        [~, SEQ{i}(:,1), ~] = fastqread(fastq_file{i,2});

        FastQData.maybe_clear_unzipped(fastq_file(i,:), ext);
        
    end
    
    Nreads = cell2mat(Nreads);
    Nreads = sum(Nreads(:,1));
    
    if (Nfastqs > 1)
        fprintf('Total reads in FASTQs: %d\n', Nreads);
    end
    
    QC = vertcat(QC{:});
    CB = vertcat(CB{:});
    UMI = vertcat(UMI{:});
    SEQ = vertcat(SEQ{:});
    
    [CB, ~, read_CB] = unique_by_freq(CB);
    [UMI, ~, read_UMI] = unique_by_freq(UMI);
    [SEQ, ~, read_SEQ] = unique_by_freq(SEQ);

end