function my_CARLIN_pipeline(SampleList,cfg_type,input_dir,output_dir,template,varargin)

    %% Note
    % we assume that the current dir is where Custom_CARLIN is
    % currently support variables: 'template', 'read_cutoff_override', 'read_cutoff_floor'
    p0 = inputParser;
    cur_dir=pwd;
    p0.addParameter('read_cutoff_override',NaN);
    p0.addParameter('read_cutoff_floor',10);
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
        sample_dir=sample_name+".trimmed.pear.assembled.fastq";

        analyze_CARLIN(char(sample_dir),char(sample_type), char(output_dir_1),'read_override_UMI_denoised',res.read_cutoff_override,'read_cutoff_UMI_denoised',res.read_cutoff_floor);

        % summary has been precomputed and saved
        output_all_from_summary(sample_name,output_dir,template,'CARLIN_dir',res.CARLIN_dir)
%         resolution=300;
%         cd(output_dir)
%         %%% plotting
%         load Summary.mat
%         close all
%         plot_highlighted_alleles(summary, length(summary.alleles)-1);
%         %     file_name="plot_highlighted_alleles.eps";
%         %     print('-depsc2','-painters',file_name);
%         file_name="highlight_alleles.png";
%         axis tight;
%         print(file_name,'-dpng',['-r' num2str(resolution)]);
%         %saveas(gcf,file_name)
% 
%         close all
%         plot_allele_frequency_CDF(summary, 'Eyeball')
%         file_name="plot_allele_frequency_CDF.png";
%         axis tight;
%         print(file_name,'-dpng',['-r' num2str(resolution)]);
%         file_name="plot_allele_frequency_CDF.eps";
%         print('-depsc2','-painters',file_name);
% 
%         close all
%         [sp, ins_freq, del_freq]= plot_indel_freq_vs_length(summary);
%         file_name="plot_indel_freq_vs_length.png";
%         axis tight;
%         print(file_name,'-dpng',['-r' num2str(resolution)]);
%         file_name="plot_indel_freq_vs_length.eps";
%         print('-depsc2','-painters',file_name);
%         save(sprintf("%s/indel_freq_vs_length.mat", output_dir), 'ins_freq', 'del_freq');
% 
%         % This is works for Tigre CARLIN data
%         close all
%         [sp, breakdown,mut_types]=plot_site_decomposition(summary, true, 'Eyeball', '# of Transcripts');
%         file_name="plot_site_decomposition.png";
%         axis tight;
%         print(file_name,'-dpng',['-r' num2str(resolution)]);
%         file_name="plot_site_decomposition.eps";
%         print('-depsc2','-painters',file_name);
%         save(sprintf("%s/site_decomposition.mat", output_dir), 'breakdown','mut_types');
% 
%         close all
%         plot_stargate.create(summary)
%         file_name="plot_stargate.png";
%         %saveas(gcf,file_name)
%         print(file_name,'-dpng',['-r' num2str(resolution)]);
%         
%         allele_freqs=summary.allele_freqs;
%         mut_list = cellfun(@(x) Mutation.identify_Cas9_events(x), summary.alleles, 'un', false); 
%         AlleleAnnotation = cellfun(@(x) arrayfun(@(i) x(i).annotate(true), [1:size(x,1)], 'un', false), mut_list, 'un', false);
%         AlleleAnnotation = cellfun(@(x) strjoin(x,','), AlleleAnnotation, 'un', false);
%         AlleleAnnotation(cellfun(@isempty, AlleleAnnotation)) = {'[]'};
%         save(sprintf("%s/allele_annotation.mat", output_dir),"allele_freqs","AlleleAnnotation");

    end
    
    % write a file to indicate that the job is done
    fileID = fopen('CARLIN_analysis_actually.done','w');
    fprintf(fileID,'Done');
    fclose(fileID);
    
    cd(cur_dir) % return to original dir
