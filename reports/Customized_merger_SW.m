function [merged_summary,sample_map,allele_breakdown_by_sample]=Customized_merger_SW(samples)

    [merged, sample_map, allele_breakdown_by_sample]=ExperimentSummary.FromMerge(samples);
    

        N.uncleaned_tags=0;
        N.cleaned_tags=0;
        N.common_tags=0;
        N.called_tags=0;
        N.eventful_tags=0;
        for jj=1:length(samples)
           N.uncleaned_tags=N.uncleaned_tags+samples(jj).N.uncleaned_tags;
           N.cleaned_tags=N.cleaned_tags+samples(jj).N.cleaned_tags;
           N.common_tags=N.common_tags+samples(jj).N.common_tags;
           N.called_tags=N.called_tags+samples(jj).N.called_tags;
           N.eventful_tags=N.eventful_tags+samples(jj).N.eventful_tags;
        end

        reads.in_fastq=0;
        reads.common_tags=0;
        reads.called_tags_total=0;
        reads.called_tags_allele=0;
        reads.eventful_tags_total=0;
        reads.eventful_tags_allele=0;
        for jj=1:length(samples)
            reads.in_fastq=reads.in_fastq+samples(jj).reads.in_fastq;
            reads.common_tags=reads.common_tags+samples(jj).reads.common_tags;
            reads.called_tags_total=reads.called_tags_total+samples(jj).reads.called_tags_total;
            reads.called_tags_allele=reads.called_tags_allele+samples(jj).reads.called_tags_allele;
            reads.eventful_tags_total=reads.eventful_tags_total+samples(jj).reads.eventful_tags_total;
            reads.eventful_tags_allele=reads.eventful_tags_allele+samples(jj).reads.eventful_tags_allele;
        end
        
        merged_summary=[];
        merged_summary.alleles=merged.alleles;
        merged_summary.allele_freqs=merged.allele_freqs;
        merged_summary.reads=reads;
        merged_summary.N=N;
        
end