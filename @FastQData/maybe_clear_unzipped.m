function maybe_clear_unzipped(fastq_file, ext)
    
    if (ischar(fastq_file))
        fastq_file = {fastq_file};
        ext = {ext};
    end
    
    for i = 1:size(fastq_file,2)
        if (strcmp(ext{i}, '.gz'))
            fprintf('Removing unzipped FASTQ file: %s\n', fastq_file{i});
            delete(fastq_file{i});            
        end
    end
end
