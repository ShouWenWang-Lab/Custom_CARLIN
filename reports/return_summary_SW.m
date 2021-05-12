function [result,annotation]=return_summary(summary,params,thresholds)

    %load Summary.mat
    %% reads
    tot_fastq_N=double(summary.reads.in_fastq);
    %common_UMI(i,j)=summary.reads.common_tags/double(summary.reads.in_fastq)*100;
    edit_read_fraction=summary.reads.eventful_tags_allele/double(summary.reads.in_fastq)*100;        

    %% UMI coverage
    x1=summary.reads.eventful_tags_total;
    x2=summary.reads.called_tags_total;
    UMI_eventful=summary.N.eventful_tags;
    UMI_called=summary.N.called_tags;
    Mean_read_per_edited_UMI=x1/UMI_eventful;
    Mean_read_per_unedited_UMI=(x2-x1)/(UMI_called-UMI_eventful);
    edit_UMI_fraction=max(round(summary.N.eventful_tags/summary.N.called_tags*100),0);



    %% Allele
    total_alleles=length(summary.alleles);
    singleton=sum(summary.allele_freqs==1);
    effective_allele_N=round(effective_alleles(summary));
    Diversity_index_all=diversity_index(summary, false); %Diversity Index (normalized by all)
    Diversity_index_edited=diversity_index(summary, true); %Diversity Index (normalized by edited)
    vector=cellfun(@(x) CARLIN_def.getInstance.N.segments-length(Mutation.find_modified_sites(x)), summary.alleles);
    CARLIN_potential_by_UMI=max(sum(vector.*summary.allele_freqs)/sum(summary.allele_freqs),0); % new form, by UMI
    CARLIN_potential_by_allel=max(mean(vector),0); % original form, by allele, not by UMI
    
    
    %% UMI
    sc = startsWith(params.Results.cfg_type, 'sc');
    
    if (sc)
        read_threshold_CB=thresholds.CB.chosen;
        read_threshold_UMI=thresholds.UMI.chosen;
        
        result=[tot_fastq_N,edit_read_fraction,read_threshold_CB,read_threshold_UMI,UMI_eventful,UMI_called,...
            Mean_read_per_edited_UMI,Mean_read_per_unedited_UMI,edit_UMI_fraction,...
            total_alleles,singleton,effective_allele_N,Diversity_index_all,Diversity_index_edited,...
            CARLIN_potential_by_UMI,CARLIN_potential_by_allel];
        annotation=["tot_fastq_N","edit_read_fraction","read_threshold_CB","read_threshold_UMI","UMI_eventful","UMI_called",...
            "Mean_read_per_edited_UMI","Mean_read_per_unedited_UMI","edit_UMI_fraction",...
            "total_alleles","singleton","effective_allele_N","Diversity_index_all","Diversity_index_edited",...
            "CARLIN_potential_by_UMI","CARLIN_potential_by_allel"];
        
    else
        read_threshold=thresholds.chosen;
        
        result=[tot_fastq_N,edit_read_fraction,read_threshold,UMI_eventful,UMI_called,...
            Mean_read_per_edited_UMI,Mean_read_per_unedited_UMI,edit_UMI_fraction,...
            total_alleles,singleton,effective_allele_N,Diversity_index_all,Diversity_index_edited,...
            CARLIN_potential_by_UMI,CARLIN_potential_by_allel];
        annotation=["tot_fastq_N","edit_read_fraction","read_threshold","UMI_eventful","UMI_called",...
            "Mean_read_per_edited_UMI","Mean_read_per_unedited_UMI","edit_UMI_fraction",...
            "total_alleles","singleton","effective_allele_N","Diversity_index_all","Diversity_index_edited",...
            "CARLIN_potential_by_UMI","CARLIN_potential_by_allel"];
    end

    
end
