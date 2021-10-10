function csv_reports(SampleList,input_dir,template)
    % the previous name: analyze_data_across_conditions_SW

    %% Note
    % we assume that the current dir is where Custom_CARLIN is
    % currently support variables: 'template', 'read_cutoff_override', 'read_cutoff_floor'
    switch_template(template)
    cur_dir=pwd;
    install_CARLIN
    
    %% start the analysis
    sample_name_array=split(SampleList,',');

    data_table={};
    annotation_array={};

    for j = 1:length(sample_name_array)

        %% prepare file and folder
        
        sample_name=sample_name_array(j);
        disp("Current sample: "+sample_name)
        
        output_dir=input_dir+"/"+sample_name;%+"_mo_v1_autoThresh3";
        mkdir(output_dir)
        cd(output_dir)
        %%% plotting
        if exist(output_dir+"/Summary.mat") == 2
            load Summary.mat
            [data_table{j},annotation_array{j}]=return_summary_SW(summary,params,thresholds);
        else
            disp("Warning: Sample"+string(sample_name)+"has not been computed")
            data_table{j}=[]; annotation_array{j}=[];

        end
    end
    
    

    %% check computed file
    max_length=0;
    max_index=0;
    for j = 1:length(sample_name_array)
       if max_length< length(data_table{j})
           max_length=length(data_table{j});
           max_index=j;
       end
    end
    
    if max_index==0
        disp("No valid files found. Abort")
    else
        
        new_folder=string(input_dir)+"/csv_report";
        mkdir(new_folder)
        cd(new_folder)
    
        annotation=annotation_array{max_index};

        for j = 1:length(sample_name_array)
            if isempty(data_table{j})
                data_table{j}=zeros(1,max_length)+NaN;
            end
        end


        matrix=cell2mat(data_table');

        %cd(input_dir)
        writematrix(matrix,'result.csv')

        fid = fopen('variable_names.txt','w');
        fid2= fopen('result.csv','a');
        for jj = 1 : length(annotation)
          fprintf( fid, annotation(jj) );
          fprintf( fid, '\n' );

          fprintf( fid2, annotation(jj) );
          fprintf( fid2, ',' );
        end
        fclose(fid);

        writematrix(matrix,'result_2.csv')


        fid = fopen('sample_names.txt','w');
        for jj = 1 : size(sample_name_array, 1)
          fprintf( fid, sample_name_array(jj) );
          fprintf( fid, '\n' );
        end
        fclose(fid);

    end
    cd(cur_dir)
    system("python generate_csv.py --path "+new_folder)
    delete(new_folder+"/result.csv")
    delete(new_folder+"/result_2.csv")
    delete(new_folder+"/sample_names.txt")
    delete(new_folder+"/variable_names.txt")
    
    
    
    
