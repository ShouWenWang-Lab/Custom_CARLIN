function my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,varargin)

    %% Note
    % we assume that the current dir is where Custom_CARLIN is
    % currently support variables: 'template', 'read_cutoff_UMI_override', 'read_cutoff_CB_override'
    p0 = inputParser;
    cur_dir=pwd;
    p0.addParameter('read_cutoff_UMI_override',NaN);
    p0.addParameter('read_cutoff_CB_override',NaN);
    p0.addParameter('CARLIN_dir',cur_dir);
    p0.parse(varargin{:});
    res=p0.Results;

    cd(res.CARLIN_dir)
    switch_template(template)
    install_CARLIN
    
    %% start the analysis
    sample_name_array=split(SampleList,',');
    sample_type_array=strings(1,length(sample_name_array));
    for j =1:length(sample_name_array)
        sample_type_array(j)=cfg_type;
    end


    for j = 1:length(sample_name_array)

        cd(input_dir)
        %    fprintf('value of a: %d\n',j)
        sample_type=sample_type_array(j);
        sample_name=sample_name_array(j);
        output_dir_1=output_dir+"/"+sample_name;
        mkdir(output_dir_1)

        if startsWith(cfg_type,'Bulk')
            sample_dir=sample_name+".trimmed.pear.assembled.fastq";

            analyze_CARLIN(char(sample_dir),char(sample_type), char(output_dir_1),'read_override_UMI_denoised',res.read_cutoff_UMI_override);
        elseif startsWith(cfg_type,'scCamellia') % we reverse R1 and R2 here
            fastq_file = {sample_name+"_R2.fastq.gz", sample_name+"_R1.fastq.gz"};
            analyze_CARLIN(fastq_file,char(sample_type), char(output_dir_1), ...
                'read_override_CB_denoised',res.read_cutoff_CB_override,'read_override_UMI_denoised',res.read_cutoff_UMI_override);
        else
            % uinsg 10X etc
            fastq_file = {sample_name+"_R1.fastq.gz", sample_name+"_R2.fastq.gz"};
            analyze_CARLIN(fastq_file,char(sample_type), char(output_dir_1), ...
                'read_override_CB_denoised',res.read_cutoff_CB_override,'read_override_UMI_denoised',res.read_cutoff_UMI_override);
        end

        % summary has been precomputed and saved
        output_all_from_summary(sample_name,output_dir,template,'CARLIN_dir',res.CARLIN_dir)

    end
    
    % write a file to indicate that the job is done
    fileID = fopen('CARLIN_analysis_actually.done','w');
    fprintf(fileID,'Done');
    fclose(fileID);
    
    cd(cur_dir) % return to original dir
