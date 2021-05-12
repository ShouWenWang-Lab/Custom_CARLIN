function [fastq_file, ext] = maybe_unzip(fastq_file)

    is_str = false;
    if (ischar(fastq_file))
        is_str = true;
        fastq_file = {fastq_file};
    end

    for i = 1:size(fastq_file,2)
        [~, file{i}, ext{i}] = fileparts(fastq_file{i});
        if (strcmp(ext{i}, '.gz'))
            unzip_path = tempname;
            fprintf('Unzipping FASTQ file to: %s/%s\n', unzip_path, file{i});
            gunzip(fastq_file{i}, unzip_path);
            fastq_file{i} = sprintf('%s/%s', unzip_path, file{i});
            assert(exist(fastq_file{i}, 'file') == 2, sprintf('Unzip unsuccessful'));
        end
        assert(endsWith(fastq_file{i}, '.fastq'), 'Unrecognized file type. Only *.fastq and *.fastq.gz supported');
    end
    
    if (is_str)
        fastq_file = fastq_file{1};
        ext = ext{1};
    end
end