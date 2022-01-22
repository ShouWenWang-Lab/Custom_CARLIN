function merge_samples(SampleList,input_dir,template,varargin)
    switch_template(template)
    p0 = inputParser;
    p0.addParameter("CARLIN_dir",".");
    p0.parse(varargin{:});
    res=p0.Results;

    cur_dir=pwd;
    cd(res.CARLIN_dir)
    install_CARLIN

    sample_name_array=split(SampleList,",");

    sample_name=sample_name_array(1);
    temp_data=load(input_dir+"/"+sample_name+"/Summary.mat");
    sample_list=[temp_data.summary];
    sample_info={};
    sample_info{1}=sample_name{1};
    if length(sample_name_array) >1
        for j = 2:length(sample_name_array)    
            sample_name=sample_name_array(j);
            temp_data=load(input_dir+"/"+sample_name+"/Summary.mat");
            sample_list(j)=temp_data.summary;
            sample_info{j}=sample_name{1};
        end
    end

    output_dir=string(input_dir)+"/merge_all";
    if ~exist(output_dir, "dir")
        % Folder does not exist so create it.
        mkdir(output_dir);
    end
    
    
   [summary,sample_map,allele_breakdown_by_sample]=Customized_merger_SW(transpose(sample_list));
    params=temp_data.params;
    thresholds=temp_data.thresholds;
    save(sprintf("%s/Summary.mat", output_dir), "summary","params","thresholds");
    sample_names=sample_name_array;

    save(sprintf("%s/allele_breakdown_by_sample.mat", output_dir),"sample_map", "allele_breakdown_by_sample","sample_names");

    allele_freqs=summary.allele_freqs;
    mut_list = cellfun(@(x) Mutation.identify_Cas9_events(x), summary.alleles, 'un', false); 
    AlleleAnnotation = cellfun(@(x) arrayfun(@(i) x(i).annotate(true), [1:size(x,1)], 'un', false), mut_list, 'un', false);
    AlleleAnnotation = cellfun(@(x) strjoin(x,','), AlleleAnnotation, 'un', false);
    AlleleAnnotation(cellfun(@isempty, AlleleAnnotation)) = {'[]'};
    save(sprintf("%s/allele_annotation.mat", output_dir),"allele_freqs","AlleleAnnotation");
    
    cd(cur_dir)