function my_CARLIN_pipeline(DataList,data_Type,data_dir,outdir_name,read_cutoff_floor,read_cutoff_override)


    if nargin<5
      read_cutoff_override = NaN;
      read_cutoff_floor=10;
    end

    if nargin<6
      read_cutoff_override = NaN;
    end

    install_CARLIN

    sample_name_array=split(DataList,',');

    ratio_by_eventful_UMI=zeros(length(sample_name_array),1);
    ratio_by_allele=zeros(length(sample_name_array),1);
    ave_insert_len_array=zeros(length(sample_name_array),1);
    ave_del_len_array=zeros(length(sample_name_array),1);
    ave_insert_del_len_ratio=zeros(length(sample_name_array),1);


    sample_type_array=strings(1,length(sample_name_array));
    for j =1:length(sample_name_array)
        sample_type_array(j)=data_Type;
    end


    for j = 1:length(sample_name_array)

        cd(data_dir)
        %    fprintf('value of a: %d\n',j)
        sample_type=sample_type_array(j);
        sample_name=sample_name_array(j);
        output_dir=outdir_name+"/"+sample_name;
        mkdir(output_dir)
        sample_dir=sample_name+".trimmed.pear.assembled.fastq";


        analyze_CARLIN(char(sample_dir),char(sample_type), char(output_dir),'read_override_UMI_denoised',read_cutoff_override,'read_cutoff_UMI_denoised',read_cutoff_floor);

        cd(data_dir)
        cd(output_dir)
        %%% plotting
        load Summary.mat
        close all
        plot_highlighted_alleles(summary, length(summary.alleles)-1);
        %     file_name="plot_highlighted_alleles.eps";
        %     print('-depsc2','-painters',file_name);
        file_name="highlight_alleles.png";
        saveas(gcf,file_name)

        close all
        plot_allele_frequency_CDF(summary, 'Eyeball')
        file_name="plot_allele_frequency_CDF.eps";
        print('-depsc2','-painters',file_name);

        close all
        plot_indel_freq_vs_length(summary)
        file_name="plot_indel_freq_vs_length.eps";
        print('-depsc2','-painters',file_name);

    % this is not working for the Tigre carlin data
    %     close all
    %     plot_site_decomposition(summary, true, 'Eyeball', '# of Transcripts')
    %     file_name="plot_site_decomposition.eps";
    %     print('-depsc2','-painters',file_name);

        close all
        plot_stargate.create(summary)
        file_name="plot_stargate.png";
        saveas(gcf,file_name)

    end