function analyze_data_across_conditions_SW(SampleList,dir_name)


    %% start the analysis
    sample_name_array=split(SampleList,',');

    ratio_by_eventful_UMI=zeros(length(sample_name_array),1);
    ratio_by_allele=zeros(length(sample_name_array),1);
    ave_insert_len_array=zeros(length(sample_name_array),1);
    ave_del_len_array=zeros(length(sample_name_array),1);

    data_table={};

    for j = 1:length(sample_name_array)

        %% prepare file and folder
        
        sample_name=sample_name_array(j);

        output_dir=dir_name+"/"+sample_name;%+"_mo_v1_autoThresh3";
        mkdir(output_dir)
        cd(output_dir)
        %%% plotting
        load Summary.mat

        [ratio_by_eventful_UMI(j),ratio_by_allele(j),ave_insert_len_array(j),ave_del_len_array(j)]=ins_del_statistics(summary);
        [data_table{j},annotation]=return_summary(summary,params,thresholds);

    end

    matrix=cell2mat(data_table');

    cd(dir_name)
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

    
