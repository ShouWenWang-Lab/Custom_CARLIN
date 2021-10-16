function make_allele_bank(SampleList,input_dir,template,varargin)
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
    bank = Bank.Create(transpose(sample_list), transpose(sample_info), output_dir);
    
    

    [summary,sample_map,allele_breakdown_by_sample]=Customized_merger_SW(transpose(sample_list));
    params=temp_data.params;
    thresholds=temp_data.thresholds;
    save(sprintf("%s/Summary.mat", output_dir), "summary","params","thresholds");
    save(sprintf("%s/allele_breakdown_by_sample.mat", output_dir),"sample_map", "allele_breakdown_by_sample","SampleList");
    %output_all_from_summary("merge_all",input_dir,template) %% need to check whether it is experimental object. Too much. 
    
    cd(cur_dir)